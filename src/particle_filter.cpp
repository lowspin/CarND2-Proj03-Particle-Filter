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

#include "particle_filter.h"

using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std_init[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 10;

  // noise generation
  default_random_engine gen;
  normal_distribution<double> N_x_init(0.0, std_init[0]);
  normal_distribution<double> N_y_init(0.0, std_init[1]);
  normal_distribution<double> N_theta_init(0.0, std_init[2]);
  Particle p;
  for (int i=0; i<num_particles; i++){
    p.id = i;
    p.x = x + N_x_init(gen);
    p.y = y + N_y_init(gen);
    p.theta = theta + N_theta_init(gen);
    p.weight = 1.0;
    particles.push_back(p);
  }
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  // noise generation for motion model
  default_random_engine gen;
  normal_distribution<double> N_x(0.0, std_pos[0]);
  normal_distribution<double> N_y(0.0, std_pos[1]);
  normal_distribution<double> N_theta(0.0, std_pos[2]);
  Particle p;
  for (int i=0; i<num_particles; i++){
    p = particles[i];
    if (fabs(yaw_rate)>1e-5){ // non-zero yaw rate
      p.x = particles[i].x + (velocity/yaw_rate)*( sin(particles[i].theta+(yaw_rate*delta_t)) - sin(particles[i].theta) ) + N_x(gen);
      p.y = particles[i].y + (velocity/yaw_rate)*( cos(particles[i].theta) - cos(particles[i].theta+(yaw_rate*delta_t)) ) + N_y(gen);
      p.theta = particles[i].theta + (yaw_rate*delta_t) + N_theta(gen);
    } else {
      p.x = particles[i].x + velocity*delta_t*cos(particles[i].theta) + N_x(gen);
      p.y = particles[i].y + velocity*delta_t*sin(particles[i].theta) + N_y(gen);
      p.theta = particles[i].theta + N_theta(gen); // no change to yaw if yaw rate is zero
    }
    particles[i] = p;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  for (int i = 0; i < observations.size(); i++) {
    observations[i].id = predicted[0].id;
    double bestdist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
    for (int j=0; j<predicted.size(); j++) {

      if (dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y)<bestdist){
        observations[i].id = predicted[j].id;
        bestdist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

  // data association
  std::vector<LandmarkObs> predicted;
  std::vector<LandmarkObs> mapped_observations;
  LandmarkObs lm;
  weights.clear();
  for (int i=0; i<num_particles; i++){

    // restrict candidates to landmarks within sensor range
    predicted.clear();
    for (int j=0; j<map_landmarks.landmark_list.size(); j++){
      lm.id = map_landmarks.landmark_list[j].id_i;
      lm.x = (double)map_landmarks.landmark_list[j].x_f;
      lm.y = (double)map_landmarks.landmark_list[j].y_f;
      if(true) { //(dist(particles[i].x, particles[i].y, lm.x, lm.y) < sensor_range){
        predicted.push_back(lm);
      }
    }

    // transform particle coordinates to map coordinates
    mapped_observations.clear();

    // Transform each particle's observations
    for(int k=0; k<observations.size(); k++) {
        lm.x = (observations[k].x * cos(particles[i].theta)) - (observations[k].y * sin(particles[i].theta)) + particles[i].x;
        lm.y = (observations[k].x * sin(particles[i].theta)) + (observations[k].y * cos(particles[i].theta)) + particles[i].y;
        mapped_observations.push_back(lm);
    }

    // find closest neighbors
    dataAssociation(predicted, mapped_observations);

    // update weights
    particles[i].weight = 1.0;
    for (int j=0; j<mapped_observations.size(); j++) {

      // locate matched landmark
      int matchid = 0;
      while (mapped_observations[j].id != predicted[matchid].id) {
        matchid++;
      }

      // calculate gaussian term
      double w = exp(
            -(pow(mapped_observations[j].x-predicted[matchid].x,2)/(2*pow(std_landmark[0],2)))
            -(pow(mapped_observations[j].y-predicted[matchid].y,2)/(2*pow(std_landmark[1],2))) )
          / (2*M_PI*std_landmark[0]*std_landmark[1]);

      // update weight with this observation
      particles[i].weight *= w;
    }

    // update global weights vector
    weights.push_back( particles[i].weight );

  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::vector<Particle> particles_new;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());

  particles_new.clear();
  for(int n=0; n<num_particles; ++n) {
      particles_new.push_back(particles[d(gen)]);
  }

  particles.clear();
  particles = particles_new;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
