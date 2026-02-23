
#include <chrono>
#include <cstdio>
#include <memory>

#include <rclcpp/logging.hpp>
#include <rclcpp/publisher.hpp>
#include <rclcpp/subscription.hpp>
#include <rclcpp/timer.hpp>
#include "rclcpp/executors.hpp"
#include "rclcpp/node.hpp"
#include "rclcpp/utilities.hpp"

#include "geometry_msgs/msg/pose_stamped.hpp"
#include "sensor_msgs/msg/point_cloud2.hpp"
#include "visualization_msgs/msg/marker_array.hpp"


#include "safe_planner/map/rog_map.hpp"
#include "pcl_conversions/pcl_conversions.h"


const std::string in_pcd_topic_name 	= "/uva/lidar/pcd"; 
const std::string in_pose_topic_name 	= "/uva/pose/stamped"; 
const std::string out_pcd_topic_name 	= "/rogmap/occupied_pcd"; 
const std::string out_bound_topic_name 	= "/rogmap/map_bound"; 

const std::string frame_id 				= "car";

class MincoTest : public rclcpp::Node {
public:
  	MincoTest() 
	: Node("map_test")
	, map_({100,100,20}, 0.1, true, 1, {0,0,0}, {})
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
	, pub_pcd_(create_publisher<sensor_msgs::msg::PointCloud2>(out_pcd_topic_name, 1))
	, pub_marker_(create_publisher<visualization_msgs::msg::MarkerArray>(out_bound_topic_name, 1))
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


			Eigen::Vector3f p,bmin,bmax;
			map_.get_local_scale(p, bmin, bmax);
			visualization_msgs::msg::MarkerArray markers;

			// 1. 无人机位置点
			{
				visualization_msgs::msg::Marker drone_marker;
				drone_marker.header.frame_id = frame_id;
				drone_marker.header.stamp = rclcpp::Clock().now();
				drone_marker.ns = "drone";
				drone_marker.id = 0;
				drone_marker.type = visualization_msgs::msg::Marker::SPHERE;
				drone_marker.action = visualization_msgs::msg::Marker::ADD;
				drone_marker.pose.position.x = p.x();
				drone_marker.pose.position.y = p.y();
				drone_marker.pose.position.z = p.z();
				drone_marker.scale.x = 0.2;
				drone_marker.scale.y = 0.2;
				drone_marker.scale.z = 0.2;
				drone_marker.color.a = 1.0;   // 不透明
				drone_marker.color.r = 1.0;   // 红色
				drone_marker.color.g = 0.0;
				drone_marker.color.b = 0.0;
				markers.markers.push_back(drone_marker);
			}

			// 2. 地图边界线框盒子
			{
				visualization_msgs::msg::Marker box_marker;
				box_marker.header.frame_id = frame_id;
				box_marker.header.stamp = rclcpp::Clock().now();
				box_marker.ns = "bound";
				box_marker.id = 1;
				box_marker.type = visualization_msgs::msg::Marker::LINE_LIST;
				box_marker.action = visualization_msgs::msg::Marker::ADD;
				box_marker.scale.x = 0.05; // 线条粗细
				box_marker.color.a = 0.5;  // 半透明
				box_marker.color.r = 0.0;
				box_marker.color.g = 1.0;  // 绿色
				box_marker.color.b = 0.0;

				// 计算边界盒子的 8 个顶点
				std::vector<Eigen::Vector3f> corners;
				corners.push_back({bmin.x(), bmin.y(), bmin.z()});
				corners.push_back({bmax.x(), bmin.y(), bmin.z()});
				corners.push_back({bmax.x(), bmax.y(), bmin.z()});
				corners.push_back({bmin.x(), bmax.y(), bmin.z()});
				corners.push_back({bmin.x(), bmin.y(), bmax.z()});
				corners.push_back({bmax.x(), bmin.y(), bmax.z()});
				corners.push_back({bmax.x(), bmax.y(), bmax.z()});
				corners.push_back({bmin.x(), bmax.y(), bmax.z()});

				// 定义盒子的 12 条边（LINE_LIST 每两个点组成一条线）
				int edges[12][2] = {
					{0,1},{1,2},{2,3},{3,0}, // 底面
					{4,5},{5,6},{6,7},{7,4}, // 顶面
					{0,4},{1,5},{2,6},{3,7}  // 竖线
				};

				for (auto &e : edges) {
					geometry_msgs::msg::Point p1, p2;
					p1.x = corners[e[0]].x(); p1.y = corners[e[0]].y(); p1.z = corners[e[0]].z();
					p2.x = corners[e[1]].x(); p2.y = corners[e[1]].y(); p2.z = corners[e[1]].z();
					box_marker.points.push_back(p1);
					box_marker.points.push_back(p2);
				}

				markers.markers.push_back(box_marker);
			}
			pub_marker_->publish(markers);
		}))
	{
		
		RCLCPP_INFO(get_logger(),"map_test_launched, run with ROGMAP");
		RCLCPP_INFO(get_logger(),"  Topics:");
		RCLCPP_INFO(get_logger(),"\t [in]  pcd   : %s", in_pcd_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [in]  pose  : %s", in_pose_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] pcd   : %s", out_pcd_topic_name.c_str());
		RCLCPP_INFO(get_logger(),"\t [out] bound : %s", out_bound_topic_name.c_str());
	}

private:
	safe_planner::map::ROGMap map_;
	rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr		sub_pcd_;
	rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr 	sub_pose_;
	
	rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr		pub_pcd_;
	rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr 	pub_marker_;

	rclcpp::TimerBase::SharedPtr timer_;
	Eigen::Vector3f robot_position_;
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
