
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

#include "safe_planner/esdf/implicit_esdf.hpp"
#include "safe_planner/planner/middle_end.hpp"

#include "geometry_msgs/msg/pose_stamped.hpp"
#include "sensor_msgs/msg/point_cloud2.hpp"
#include "nav_msgs/msg/path.hpp"
#include "std_msgs/msg/float32_multi_array.hpp"

#include <pcl_conversions/pcl_conversions.h>

#include "safe_planner/map/rog_map.hpp"
#include "safe_planner/planner/front_end.hpp"
#include "safe_planner/trajectory/multi_poly.hpp"
#include "visualization_msgs/msg/marker.hpp"
#include "visualization_msgs/msg/marker_array.hpp"



const std::string in_pcd_topic_name 	= "/uva/lidar/pcd"; 
const std::string in_pose_topic_name 	= "/uva/pose/stamped";
const std::string in_target_topic_name 	= "/sim/target"; 

const std::string out_trajectory_topic_name 		= "/middle_end/trajectory"; 
const std::string out_pcd_topic_name 				= "/rogmap/occupied_pcd"; 
const std::string out_vel_topic_name 				= "/middle_end/velocity"; 
const std::string out_acc_topic_name 				= "/middle_end/acceleration";
const std::string out_jerk_topic_name 				= "/middle_end/jerk"; 
const std::string out_snap_topic_name 				= "/middle_end/snap";
const std::string out_time_topic_name 				= "/middle_end/time"; 
const std::string out_path_topic_name 				= "/front_end/path"; 
const std::string out_control_points_topic_name 	= "/middle_end/control_points";

const std::string frame_id 				= "car";

class MincoTest : public rclcpp::Node {
public:
  	MincoTest() 
	: Node("minco_test")
	, map_({200,200,20}, 0.1, true, 1, {0,0,0}, {})
	, front_end_(map_,{.algo = safe_planner::planner::front_end::Config::Algo::RRT})
	, esdf_(map_, 0.1)
	, middle_end_(map_,esdf_)
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
	, pub_time_(create_publisher<std_msgs::msg::Float32MultiArray>(out_time_topic_name, 1))
	, pub_vel_(create_publisher<std_msgs::msg::Float32MultiArray>(out_vel_topic_name, 1))
	, pub_acc_(create_publisher<std_msgs::msg::Float32MultiArray>(out_acc_topic_name, 1))
	, pub_jerk_(create_publisher<std_msgs::msg::Float32MultiArray>(out_jerk_topic_name, 1))
	, pub_snap_(create_publisher<std_msgs::msg::Float32MultiArray>(out_snap_topic_name, 1))
	, pub_control_points_(create_publisher<visualization_msgs::msg::MarkerArray>(out_control_points_topic_name, 1))
	, tt_(1)
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

			static std::vector<Eigen::Vector3f> points{
				{0,0,0},
				{1,0,0},
				{1,1,0},
				{1,1,1}
			};
			static std::vector<Eigen::Vector3d> pointsd{};
			static std::vector<Eigen::Vector3d> control_points;
			static std::vector<Eigen::Vector3d> control_gradients;

			safe_planner::planner::trajectory::UniformBSpline trajectory;
			std::string info;
			front_end_.search(robot_position_, robot_target_, points);
			pointsd.clear();
			for(const auto& p : points){
				pointsd.emplace_back(p.cast<double>());
			}
			tt_.Begin();
			middle_end_.optimize(
				pointsd, 
				Eigen::Matrix<double, 3, 2>::Zero(), 
				Eigen::Matrix<double, 3, 2>::Zero(), 
				trajectory);
			tt_.End();
			if(!tt_.Log(info)){
				return;
			}
			RCLCPP_INFO(get_logger(),"[MiddleEnd Run] %s",info.c_str());
			// middle_end_.test_optimizer(trajectory, control_gradients);
			std_msgs::msg::Float32MultiArray time_msg;
			std_msgs::msg::Float32MultiArray vel_msg;
			std_msgs::msg::Float32MultiArray acc_msg;
			std_msgs::msg::Float32MultiArray jerk_msg;
			std_msgs::msg::Float32MultiArray snap_msg;
			nav_msgs::msg::Path traj_msg;
			traj_msg.header.frame_id = frame_id;
			safe_planner::planner::ITrajectory::seconds t1,t2,err{0.01f};
			trajectory.get_time_range(t1,t2);
			for(auto t = t1; t < t2; t += err){
				Eigen::Vector3d p = trajectory.pos(t);
				time_msg.data.push_back(t.count());
				vel_msg.data.push_back(trajectory.vel(t));
				acc_msg.data.push_back(trajectory.acc(t));
				jerk_msg.data.push_back(trajectory.jerk(t));
				snap_msg.data.push_back(trajectory.snap(t));
				geometry_msgs::msg::PoseStamped pose;
				pose.header.frame_id = frame_id;
				pose.pose.position.x = p.x();
				pose.pose.position.y = p.y();
				pose.pose.position.z = p.z();
				traj_msg.poses.push_back(pose);
			}
			time_msg.layout.dim.push_back(std_msgs::msg::MultiArrayDimension());
			time_msg.layout.dim[0].size = time_msg.data.size();
			time_msg.layout.dim[0].stride = time_msg.data.size();
			time_msg.layout.dim[0].label = "time";
			vel_msg.layout.dim.push_back(std_msgs::msg::MultiArrayDimension());
			vel_msg.layout.dim[0].size = vel_msg.data.size();
			vel_msg.layout.dim[0].stride = vel_msg.data.size();
			vel_msg.layout.dim[0].label = "velocity";
			acc_msg.layout.dim.push_back(std_msgs::msg::MultiArrayDimension());
			acc_msg.layout.dim[0].size = acc_msg.data.size();
			acc_msg.layout.dim[0].stride = acc_msg.data.size();
			acc_msg.layout.dim[0].label = "acceleration";
			jerk_msg.layout.dim.push_back(std_msgs::msg::MultiArrayDimension());
			jerk_msg.layout.dim[0].size = acc_msg.data.size();
			jerk_msg.layout.dim[0].stride = acc_msg.data.size();
			jerk_msg.layout.dim[0].label = "jerk";
			snap_msg.layout.dim.push_back(std_msgs::msg::MultiArrayDimension());
			snap_msg.layout.dim[0].size = acc_msg.data.size();
			snap_msg.layout.dim[0].stride = acc_msg.data.size();
			snap_msg.layout.dim[0].label = "snap";
			pub_traj_->publish(traj_msg);
			pub_time_->publish(time_msg);
			pub_vel_->publish(vel_msg);
			pub_acc_->publish(acc_msg);
			pub_jerk_->publish(jerk_msg);
			pub_snap_->publish(snap_msg);

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
			// trajectory.control_point(control_points);
			// static visualization_msgs::msg::MarkerArray control_points_msg{};
			// for(size_t i = 0;i < control_points_msg.markers.size();++i)
			// 	control_points_msg.markers[i].action = visualization_msgs::msg::Marker::DELETEALL;
			// pub_control_points_->publish(control_points_msg);
			// control_points_msg.markers.clear();
			// for(size_t i = 0; i < control_points.size(); ++i){
			// 	visualization_msgs::msg::Marker pose;
			// 	const auto p = control_points[i];
				// const auto g = control_gradients[i];
				// pose.header.frame_id = frame_id;
				// pose.pose.position.x = p.x();
				// pose.pose.position.y = p.y();
				// pose.pose.position.z = p.z();
				// Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitX(), g.normalized());
				// pose.pose.orientation.x = q.x();
				// pose.pose.orientation.y = q.y();
				// pose.pose.orientation.z = q.z();
				// pose.pose.orientation.w = q.w();
			// 	pose.id = i;
			// 	pose.type = visualization_msgs::msg::Marker::ARROW;
			// 	// pose.scale.x = g.norm();
			// 	pose.scale.y = 0.1;
			// 	pose.scale.z = 0.1;
			// 	pose.lifetime.sec = 1;
			// 	control_points_msg.markers.emplace_back(pose);
			// }
			// pub_control_points_->publish(control_points_msg);
		}))
	{
		
		RCLCPP_INFO(get_logger(),"rrt_star, run with ROGMAP");
		RCLCPP_INFO(get_logger(),"  Topics:");
		RCLCPP_INFO(get_logger(),"\t [in]  pcd   : %s", in_pcd_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [in]  pose  : %s", in_pose_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [in]  target: %s", in_target_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] time  : %s", out_time_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] vel   : %s", out_vel_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] acc   : %s", out_acc_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] jerk  : %s", out_jerk_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] snap  : %s", out_snap_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] pcd   : %s", out_pcd_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] path  : %s", out_path_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] traj  : %s", out_trajectory_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] control points  : %s", out_control_points_topic_name.c_str());
	}

private:
	std::vector<Eigen::Vector3d> path_{};

	safe_planner::map::ROGMap map_;
	safe_planner::planner::FrontEnd<safe_planner::map::ROGMap> front_end_;
	safe_planner::esdf::ImplicitESDF<safe_planner::map::ROGMap> esdf_;
	safe_planner::planner::MiddleEnd<safe_planner::map::ROGMap> middle_end_;

	rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr		sub_pcd_;
	rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr 	sub_pose_;
	rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr 	sub_target_;

	
	rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr		pub_pcd_;
	rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr				pub_path_;
	rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr				pub_traj_;
	rclcpp::Publisher<std_msgs::msg::Float32MultiArray>::SharedPtr	pub_time_;
	rclcpp::Publisher<std_msgs::msg::Float32MultiArray>::SharedPtr	pub_vel_;
	rclcpp::Publisher<std_msgs::msg::Float32MultiArray>::SharedPtr	pub_acc_;
	rclcpp::Publisher<std_msgs::msg::Float32MultiArray>::SharedPtr	pub_jerk_;
	rclcpp::Publisher<std_msgs::msg::Float32MultiArray>::SharedPtr	pub_snap_;
	rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_control_points_;

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
