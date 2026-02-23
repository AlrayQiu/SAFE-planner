
#include <chrono>
#include <cstdio>
#include <memory>
#include <vector>

#include <rclcpp/logging.hpp>
#include <rclcpp/publisher.hpp>
#include <rclcpp/subscription.hpp>
#include <rclcpp/timer.hpp>
#include "rclcpp/executors.hpp"
#include "rclcpp/node.hpp"
#include "rclcpp/utilities.hpp"

#include "minco_test/time_test.hpp"

#include "safe_planner/planner/middle_end.hpp"
#include "geometry_msgs/msg/pose_stamped.hpp"
#include "sensor_msgs/msg/point_cloud2.hpp"
#include "nav_msgs/msg/path.hpp"

#include <pcl_conversions/pcl_conversions.h>

#include "safe_planner/map/rog_map.hpp"
#include "safe_planner/planner/front_end.hpp"
#include "safe_planner/trajectory/minco.hpp"



const std::string in_pcd_topic_name 	= "/uva/lidar/pcd"; 
const std::string in_pose_topic_name 	= "/uva/pose/stamped";
const std::string in_target_topic_name 	= "/sim/target"; 

const std::string out_trajectory_topic_name 		= "/middle_end/trajectory"; 
const std::string out_pcd_topic_name 				= "/rogmap/occupied_pcd"; 
const std::string out_path_topic_name 				= "/front_end/path"; 
const std::string out_control_points_topic_name 	= "/middle_end/control_points";

const std::string frame_id 				= "car";

class MincoTest : public rclcpp::Node {
public:
  	MincoTest() 
	: Node("rrt_star_test")
	, map_({100,100,20}, 0.1, true, 1, {0,0,0}, {})
	, front_end_(map_,{.algo = safe_planner::planner::front_end::Config::Algo::RRT})
	, middle_end_(map_)
	, sub_pcd_(
		create_subscription<sensor_msgs::msg::PointCloud2>(
		in_pcd_topic_name, 10, 
		[this](const sensor_msgs::msg::PointCloud2::SharedPtr msg){
			pcd_.clear();
			pcl::fromROSMsg(*msg, pcd_);
		}))
	, sub_pose_(
		create_subscription<geometry_msgs::msg::PoseStamped>(
		in_pose_topic_name, 10, 
		[this](const geometry_msgs::msg::PoseStamped::SharedPtr msg){
			robot_position_(0) = msg->pose.position.x;
			robot_position_(1) = msg->pose.position.y;
			robot_position_(2) = msg->pose.position.z;
		}))
	, sub_target_(
		create_subscription<geometry_msgs::msg::PoseStamped>(
		in_target_topic_name, 10,
		[this](const geometry_msgs::msg::PoseStamped::SharedPtr msg){
			robot_target_(0) = msg->pose.position.x;
			robot_target_(1) = msg->pose.position.y;
			robot_target_(2) = msg->pose.position.z;
		}))
	, pub_pcd_(create_publisher<sensor_msgs::msg::PointCloud2>(out_pcd_topic_name, 1))
	, pub_path_(create_publisher<nav_msgs::msg::Path>(out_path_topic_name, 1))
	, pub_traj_(create_publisher<nav_msgs::msg::Path>(out_trajectory_topic_name, 1))
	, pub_control_points_(create_publisher<nav_msgs::msg::Path>(out_control_points_topic_name, 1))
	, tt_()
	, timer_(
		create_timer(std::chrono::milliseconds(100), 
		[this](){
			map_.update(pcd_, robot_position_);
			
			pcl::PointCloud<pcl::PointXYZI> occupied_pcd{};
			sensor_msgs::msg::PointCloud2 	msg{};
			map_.get_occupied_points(occupied_pcd);
			pcl::toROSMsg(occupied_pcd, msg);
			msg.header.frame_id = frame_id;
			pub_pcd_->publish(msg);
			occupied_pcd.clear();

			static std::vector<Eigen::Vector3f> points;
			static std::vector<Eigen::Vector3f> control_points;

			safe_planner::planner::trajectory::Minco trajectory;
			std::string info;
			front_end_.search(robot_position_, robot_target_, points);
			tt_.Begin();
			middle_end_.optimize(points, Eigen::Matrix<float, 2, 3>::Zero(), Eigen::Matrix<float, 2, 3>::Zero(), trajectory);
			tt_.End();
			if(tt_.Log(info)){
				RCLCPP_INFO(get_logger(),"[RRTStar Run] %s",info.c_str());
			}
			nav_msgs::msg::Path traj_msg;
			traj_msg.header.frame_id = frame_id;
			safe_planner::planner::ITrajectory::seconds t1,t2,err{0.1f};
			trajectory.get_time_range(t1,t2);
			for(auto t = t1; t < t2; t += err){
				Eigen::Vector3f p = trajectory.pos(safe_planner::planner::trajectory::Minco::seconds{t});
				geometry_msgs::msg::PoseStamped pose;
				pose.header.frame_id = frame_id;
				pose.pose.position.x = p.x();
				pose.pose.position.y = p.y();
				pose.pose.position.z = p.z();
				traj_msg.poses.push_back(pose);
			}
			pub_traj_->publish(traj_msg);

			nav_msgs::msg::Path path_msg;
			path_msg.header.frame_id = frame_id;
			for(const auto& p : points){
				geometry_msgs::msg::PoseStamped pose;
				pose.header.frame_id = frame_id;
				pose.pose.position.x = p.x();
				pose.pose.position.y = p.y();
				pose.pose.position.z = p.z();
				path_msg.poses.push_back(pose);
			}
			pub_path_->publish(path_msg);
			trajectory.control_point(control_points);
			nav_msgs::msg::Path control_points_msg;
			control_points_msg.header.frame_id = frame_id;
			for(const auto& p : control_points){
				geometry_msgs::msg::PoseStamped pose;
				pose.header.frame_id = frame_id;
				pose.pose.position.x = p.x();
				pose.pose.position.y = p.y();
				pose.pose.position.z = p.z();
				control_points_msg.poses.push_back(pose);
			}
			pub_control_points_->publish(control_points_msg);
		}))
	{
		
		RCLCPP_INFO(get_logger(),"rrt_star, run with ROGMAP");
		RCLCPP_INFO(get_logger(),"  Topics:");
		RCLCPP_INFO(get_logger(),"\t [in]  pcd   : %s", in_pcd_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [in]  pose  : %s", in_pose_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [in]  target: %s", in_target_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] pcd   : %s", out_pcd_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] path  : %s", out_path_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] traj  : %s", out_trajectory_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] control points  : %s", out_control_points_topic_name.c_str());
	}

private:
	std::vector<Eigen::Vector3f> path_{};

	safe_planner::map::ROGMap map_;
	safe_planner::planner::FrontEnd<safe_planner::map::ROGMap> front_end_;
	safe_planner::planner::MiddleEnd<safe_planner::map::ROGMap, safe_planner::planner::trajectory::Minco> middle_end_;

	rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr		sub_pcd_;
	rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr 	sub_pose_;
	rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr 	sub_target_;

	
	rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr		pub_pcd_;
	rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr				pub_path_;
	rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr				pub_traj_;
	rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr				pub_control_points_;

	TimeTest tt_;
	rclcpp::TimerBase::SharedPtr timer_;
	Eigen::Vector3f robot_position_;
	Eigen::Vector3f robot_target_;
	pcl::PointCloud<pcl::PointXYZINormal> pcd_;
};

int main(int argc, char **argv) {
  	(void)argc;
  	(void)argv;
  	rclcpp::init(argc, argv);
  	auto node = std::make_shared<MincoTest>();
  	rclcpp::spin(node);
  	rclcpp::shutdown();
}
