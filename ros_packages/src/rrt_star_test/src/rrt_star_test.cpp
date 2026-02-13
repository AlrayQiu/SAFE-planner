
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

#include "geometry_msgs/msg/pose_stamped.hpp"
#include "rrt_star_test/time_test.hpp"
#include "sensor_msgs/msg/point_cloud2.hpp"
#include "nav_msgs/msg/path.hpp"
#include "pcl_conversions/pcl_conversions.h"


#include "safe_planner/map/rog_map.hpp"
#include "safe_planner/planner/front_end.hpp"



const std::string in_pcd_topic_name 	= "/uva/lidar/pcd"; 
const std::string in_pose_topic_name 	= "/uva/pose/stamped"; 
const std::string in_target_topic_name 	= "/sim/target"; 
const std::string out_pcd_topic_name 	= "/rogmap/occupied_pcd"; 
const std::string out_path_topic_name 	= "/front_end/path"; 

const std::string frame_id 				= "car";

class RRTStarTest : public rclcpp::Node {
public:
  	RRTStarTest() 
	: Node("rrt_star_test")
	, map_({100,100,20}, 0.1, true, 1, {0,0,0}, {})
	, front_end_(map_,{.algo = safe_planner::planner::front_end::Config::Algo::RRT})
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
			tt_.Begin();
			front_end_.search(robot_position_, robot_target_, path_);
			tt_.End();
			std::string info;
			if(tt_.Log(info)){
				RCLCPP_INFO(get_logger(),"[RRTStar Run] %s",info.c_str());
			}
			nav_msgs::msg::Path path{};
			path.header.frame_id = frame_id;

			for (const auto& point : path_) {
				geometry_msgs::msg::PoseStamped pose;
				pose.header = path.header;
				pose.pose.position.x = point.x();
				pose.pose.position.y = point.y();
				pose.pose.position.z = point.z();

				pose.pose.orientation.w = 1.0;
				pose.pose.orientation.x = 0.0;
				pose.pose.orientation.y = 0.0;
				pose.pose.orientation.z = 0.0;

				path.poses.emplace_back(std::move(pose));
			}
			pub_path_->publish(path);
		}))
	{
		
		RCLCPP_INFO(get_logger(),"rrt_star, run with ROGMAP");
		RCLCPP_INFO(get_logger(),"  Topics:");
		RCLCPP_INFO(get_logger(),"\t [in]  pcd   : %s", in_pcd_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [in]  pose  : %s", in_pose_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [in]  tar   : %s", in_target_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] pcd   : %s", out_pcd_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] path  : %s", out_path_topic_name.c_str());
	}

private:
	std::vector<Eigen::Vector3f> path_{};

	safe_planner::map::ROGMap map_;
	safe_planner::planner::FrontEnd<safe_planner::map::ROGMap> front_end_;

	rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr		sub_pcd_;
	rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr 	sub_pose_;
	rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr 	sub_target_;
	
	rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr		pub_pcd_;
	rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr				pub_path_;

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
  	auto node = std::make_shared<RRTStarTest>();
  	rclcpp::spin(node);
  	rclcpp::shutdown();
}
