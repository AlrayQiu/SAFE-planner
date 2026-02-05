from launch import LaunchDescription
from launch_ros.actions import Node

def generate_launch_description():

    # Unity TCP Endpoint（ROS–Unity 通信）
    unity_endpoint = Node(
        package="ros_tcp_endpoint",
        executable="default_server_endpoint",
        name="unity_tcp_endpoint",
        output="screen",
        parameters=[{
            "ROS_IP": "0.0.0.0",     # 监听所有地址
            "ROS_TCP_PORT": 40425    # Unity 默认端口
        }]
    )

    # Foxglove Bridge
    foxglove = Node(
        package="foxglove_bridge",
        executable="foxglove_bridge",
        name="foxglove_bridge",
        output="screen",
        parameters=[{
            "port": 8765,
            "address": "0.0.0.0",
            "send_buffer_limit" : 104857600,
        }]
    )

    return LaunchDescription([
        unity_endpoint,
        foxglove
    ])
