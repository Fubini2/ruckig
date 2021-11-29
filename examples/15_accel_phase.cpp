#include <iostream>

#include <ruckig/ruckig.hpp>


using namespace ruckig;

int main() {
    // Create input parameters
    InputParameter<3> input;
    input.current_position = { 30.0, -50.0, 40.0 };
    input.current_velocity = { 10.0, -10.0, 10.0 };
    input.current_acceleration = { 10.0, -10.0, -10.0 };

    input.target_position = { -30.0, 50.0, -40.0 };
    input.target_velocity = { -10.0, -10.0, 10.0 };
    input.target_acceleration = { 10.0, 10.0, -10.0 };

    input.max_velocity = { 5.0, 5.0, 5.0 };
    input.max_acceleration = { 5.0, 5.0, 5.0 };
    input.max_jerk = { 100.0, 100.0, 100.0 };

    // We don't need to pass the control rate (cycle time) when using only offline features
    Ruckig<3> otg;
    Trajectory<3> trajectory;

    // Calculate the trajectory in an offline manner (outside of the control loop)
    Result result = otg.calculate(input, trajectory);
    if (result == Result::ErrorInvalidInput) {
        std::cout << "Invalid input!" << std::endl;
        return -1;
    }

    // Get duration of the trajectory
    std::cout << "Trajectory duration: " << trajectory.get_duration() << " [s]." << std::endl;

    double new_time{ 30.0 };

    // Then, we can calculate the kinematic state at a given time
    std::array<double, 3> new_position, new_velocity, new_acceleration;
    trajectory.at_time(new_time, new_position, new_velocity, new_acceleration);

    std::cout << "Position at time " << new_time << " [s]: " << new_position[0] << ", " << new_position[1] << ", " << new_position[2] << std::endl;

    // Get some info about the position extrema of the trajectory
    std::array<Bound, 3> position_extrema = trajectory.get_position_extrema();
    std::cout << "Position extremas for DoF 1 are " << position_extrema[0].min << " (min) to " << position_extrema[0].max << " (max)" << std::endl;
}
