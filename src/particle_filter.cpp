/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		Particle particle;
		particle.id = i;
		particle.x = dist_x(random_gen); 
		particle.y = dist_y(random_gen);
		particle.theta = dist_theta(random_gen);
		particle.weight = 1.0;
		particles.push_back(particle);
		weights.push_back(1.0);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	normal_distribution<double> dist_x(0.0, std_pos[0]); 
	normal_distribution<double> dist_y(0.0, std_pos[1]);
	normal_distribution<double> dist_theta(0.0, std_pos[2]);
	for (int i = 0; i < num_particles; ++i) {
		if (yaw_rate == 0.0) {
		particles[i].x += velocity * delta_t * cos(particles[i].theta);
		particles[i].y += velocity * delta_t * sin(particles[i].theta);
		} else {
		particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
		particles[i].y += velocity / yaw_rate * (cos(particles[i].theta)-cos(particles[i].theta + yaw_rate * delta_t));
		particles[i].theta += yaw_rate * delta_t;
		}
		// add noise
		particles[i].x += dist_x(random_gen);
		particles[i].y += dist_y(random_gen);
		particles[i].theta += dist_theta(random_gen);
	}	

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	//const double maxErr = pow(sensor_range, 2);
	const std::vector<Map::single_landmark_s> &landmarks = map_landmarks.landmark_list;
	for (int i = 0; i < num_particles; ++i) {
		Particle &particle = particles[i];

		double R[2][2];
		R[0][0] = cos(particle.theta); R[0][1] = -sin(particle.theta);
		R[1][0] = sin(particle.theta); R[1][1] = cos(particle.theta);
		double p[2];
		p[0] = particle.x; p[1] = particle.y;
		particle.sense_x.clear(); particle.sense_y.clear();
		for (const LandmarkObs &observation: observations) {
			particle.sense_x.push_back(R[0][0]*observation.x + R[0][1]*observation.y + p[0]);
			particle.sense_y.push_back(R[1][0]*observation.x + R[1][1]*observation.y + p[1]);
		}
		particle.associations.clear();
		for (int i = 0; i < observations.size(); ++i) {
				double minDist = numeric_limits<double>::max();
				int minIdx = 0;
				for (int j = 0; j < landmarks.size(); ++j) {
					double dist = (
					pow(particle.sense_x[i] - landmarks[j].x_f, 2) + pow(particle.sense_y[i] - landmarks[j].y_f, 2)
				);
					if (dist < minDist) {
						minDist = dist;
						minIdx = j;
					}
				}
				particle.associations.push_back(minIdx + 1);
		}

		for (int j = 0; j < particle.associations.size(); ++j) {
			const double xSqrErr = pow(particle.sense_x[j] - landmarks[particle.associations[j] -1].x_f, 2);
			const double ySqrErr = pow(particle.sense_y[j] - landmarks[particle.associations[j] -1].y_f, 2);
/*
			if (xSqrErr > maxErr || ySqrErr > maxErr ||xSqrErr + ySqrErr > maxErr) {
				weights[i] = 0.0;
				break;
			} else {
*/
			double weight_tem = 0.5/(M_PI*std_landmark[0]*std_landmark[1])*
								exp(-0.5*(xSqrErr/pow(std_landmark[0], 2) + ySqrErr/pow(std_landmark[1], 2)));
			weights[i] *= weight_tem;
			//}
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: Y	const std::vector<Map::single_landmark_s> &landmarks = map_landmarks.landmark_list;ou may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	discrete_distribution<> newParticleIdx(weights.begin(), weights.end());

	vector<Particle> newParticles;
	vector<double> newWeights;
	for (int i = 0; i < num_particles; ++i) {
		newParticles.push_back(particles[newParticleIdx(random_gen)]);
		newWeights.push_back(1.0);
	}

	particles = newParticles;
	weights = newWeights;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
