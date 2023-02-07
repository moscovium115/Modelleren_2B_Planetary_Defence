//
// Created by anis on 13-12-22.
//
//
// Created by anis on 30-11-22.
//
#include "N_body.h"
#include "cmath"
#include <fstream>
#include <iostream>
#include "iomanip"


#include <vector>
#include <bits/stdc++.h>

using namespace std;
void print(vector<vector<double>> &listOfVectors)
{
    for (auto vect : listOfVectors) {
        // Each element of the list is
        // a vector itself
        vector<double> currentVector = vect;

        cout << "[ ";

        // Printing vector contents
        for (auto element : currentVector)
            cout << element << ' ';

        cout << ']';
        cout << '\n';
    }
}

vector<vector<double> > perturbations(vector<double>& start_pos, double r_unc, double d_r, double d_theta, double d_phi) {
    vector<vector<double> > listOfVectors;
    double x = start_pos[0];
    double y = start_pos[1];
    double z = start_pos[2];

    for(double r = d_r; r <= r_unc; r += d_r ){
        for(double theta = 0; theta <= M_PI; theta += d_theta) {
            for(double phi = 0; phi <= 2*M_PI; phi += d_phi) {
                double dx = r*sin(theta)*cos(phi);
                double dy = r*sin(theta)*sin(phi);
                double dz = r*cos(theta);

                vector<double> vec_car;

                vec_car.push_back(x+dx);
                vec_car.push_back(y+dy);
                vec_car.push_back(z+dz);

                listOfVectors.push_back(vec_car);

            }

        }


    }
    return listOfVectors;
}





//Print vector
void print(std::vector<float> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        std::cout << input.at(i) << ' ';
    }
    std::cout << ' ' << std::endl;
}

void print_double(std::vector<double> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        std::cout << input.at(i) << ' ';
    }
    std::cout << ' ' << std::endl;
}

void Universe::use_initial_conditions_func(double &t){
    ifstream inFile;
    inFile.open("final_state.txt");
    if (!inFile) {
        cerr << "Unable to open file datafile.txt";
        exit(1);   // call system to stop
    }
    double x;
    std::vector<double> inital_values;
    while (inFile >> x) {
//        my_str.erase(remove(my_str.begin(), my_str.end(), 'A'), my_str.end());
        inital_values.push_back(x);
    }
    inFile.close();
    for (int i = 0; i < planets.size(); i++) {
        auto &planet = planets[i];
        planet.x_pos = inital_values[9 * i];
        planet.y_pos = inital_values[9 * i + 1];
        planet.z_pos = inital_values[9 * i + 2];
        planet.x_vel = inital_values[9 * i+3];
        planet.y_vel = inital_values[9 * i + 4];
        planet.z_vel = inital_values[9 * i + 5];
        planet.x_acc = inital_values[9 * i+6];
        planet.y_acc = inital_values[9 * i + 7];
        planet.z_acc = inital_values[9 * i + 8];
//        std::cout<<"x pos planeet: "<<planet.x_pos<<std::endl;
    }
    t = inital_values[6 * planets.size()];

}

void Universe::energies(std::vector<double> &total_impulse,std::vector<double> &ang_momentum, double &sum_kin_energy){
    total_impulse={0,0,0};
    sum_kin_energy=0;
    ang_momentum={0,0,0};
    for (int i = 0; i < planets.size(); i++) {
        auto &planet = planets[i];
        sum_kin_energy += 0.5 * planet.mass *(planet.x_vel * planet.x_vel + planet.y_vel * planet.y_vel + planet.z_vel * planet.z_vel);
        total_impulse[0] += planet.mass * planet.x_vel;
        total_impulse[1] += planet.mass * planet.y_vel;
        total_impulse[2] += planet.mass * planet.z_vel;
        ang_momentum[0] +=planet.mass * (planet.y_vel * planet.z_pos - planet.z_pos * planet.y_vel);
        ang_momentum[1] +=planet.mass * (planet.z_vel * planet.x_pos - planet.z_pos * planet.x_vel);
        ang_momentum[2] +=planet.mass * (planet.y_vel * planet.x_pos - planet.y_pos * planet.x_vel);

    }
    std::cout <<  " Total Impulse: " << std::endl;
    print_double(total_impulse);
    std::cout <<  " angular momentum: " << std::endl;
    print_double(ang_momentum);
}

void Universe::update_dt(double dt, bool dt_changed){
    auto &Earth=planets[1];
    auto &asteroid=planets[2];
    double distance_asteroide=sqrt((Earth.x_pos-asteroid.x_pos)*(Earth.x_pos-asteroid.x_pos)+(Earth.y_pos-asteroid.y_pos)*(Earth.y_pos-asteroid.y_pos)+(Earth.z_pos-asteroid.z_pos)*(Earth.z_pos-asteroid.z_pos));

    if (distance_asteroide<1e10 and not dt_changed){
        dt=dt/10;
        dt_changed=true;
    }
    if (dt_changed and distance_asteroide>1e10){
        dt=dt*10;
        dt_changed=false;


    }
}

void Universe::detect_collision(){
    auto &Earth=planets[1];
    auto &asteroid=planets[2];
    double distance_asteroide=sqrt((Earth.x_pos-asteroid.x_pos)*(Earth.x_pos-asteroid.x_pos)+(Earth.y_pos-asteroid.y_pos)*(Earth.y_pos-asteroid.y_pos)+(Earth.z_pos-asteroid.z_pos)*(Earth.z_pos-asteroid.z_pos));
    if (distance_asteroide<1e11){
//    std::cout<<"distance earth asteroide " << distance_asteroide<<std::endl;
    }
    if (distance_asteroide<Earth.radius){
        std::cout<<"Earth has been hit!!!!!"<<std::endl;
        std::cout<<"distance earth asteroide " << distance_asteroide<<std::endl;
//        std::cout<<"aard x pos "<<Earth.x_pos<<std::endl;
//        std::cout<<"aard y pos "<<Earth.y_pos<<std::endl;
//        std::cout<<"aard z pos "<<Earth.z_pos<<std::endl;
//        std::cout<<"aard x vel "<<Earth.x_vel<<std::endl;
//        std::cout<<"aard y vel "<<Earth.y_vel<<std::endl;
//        std::cout<<"aard z vel "<<Earth.z_vel<<std::endl;


    }
}

void Universe::detect_collision_2(int &botsing_aantal){
    auto &Earth=planets[1];
    auto &asteroid=planets[2];
    double distance_asteroide=sqrt((Earth.x_pos-asteroid.x_pos)*(Earth.x_pos-asteroid.x_pos)+(Earth.y_pos-asteroid.y_pos)*(Earth.y_pos-asteroid.y_pos)+(Earth.z_pos-asteroid.z_pos)*(Earth.z_pos-asteroid.z_pos));
    if (distance_asteroide<1e11){
//    std::cout<<"distance earth asteroide " << distance_asteroide<<std::endl;
    }
    if (distance_asteroide<Earth.radius){
        std::cout<<"Earth has been hit!!!!!"<<std::endl;
        std::cout<<"distance earth asteroide " << distance_asteroide<<std::endl;
        botsing_aantal+=1;
    }
}

void Universe::inelastic_collision(int i) {
    if (i == 10) {
        std::cout << "dart!!!!" << std::endl;
        double mass_dart = 1e20;
        double speed_dart_x = -1e5;
        double speed_dart_y = 0;
        double speed_dart_z = 0;
        auto &asteroid = planets[4];
        asteroid.x_vel = (asteroid.x_vel * asteroid.mass + speed_dart_x * mass_dart) / (asteroid.mass + mass_dart);
        asteroid.y_vel = (asteroid.y_vel * asteroid.mass + speed_dart_y * mass_dart) / (asteroid.mass + mass_dart);
        asteroid.z_vel = (asteroid.z_vel * asteroid.mass + speed_dart_z * mass_dart) / (asteroid.mass + mass_dart);
        asteroid.mass += mass_dart;
    }
}


void Universe::save_state(double &t) {
    std::cout << "saving final state" << std::endl;
    ofstream myfile;
    myfile.open("final_state.txt");
    for (int astastsa = 0; astastsa < planets.size(); astastsa++) {
        auto &planet = planets[astastsa];
        myfile << std::setw(30) << std::setprecision(100) << std::fixed << planet.x_pos << "\n " << planet.y_pos << "\n" << planet.z_pos << "\n" << planet.x_vel << "\n" << planet.y_vel
               << "\n" << planet.z_vel << "\n"<< planet.x_acc<< "\n" << planet.y_acc << "\n" << planet.z_acc << "\n";
        std::cout<<"Saving positions... "<<"\n"<<setprecision(19)<<planet.x_pos<<" "<<planet.y_pos<<" "<<planet.z_pos<<std::endl;
        std::cout<<"Saving velocities... "<<"\n"<<setprecision(19)<<planet.x_vel<<" "<<planet.y_vel<<" "<<planet.z_vel<<std::endl;

    }
    myfile << t<<"\n";
    myfile.close();
}



//geef i ook als arugment check wanneer i=0, dan doe je die beginvoorwaarde voor v0.5
//void Universe::step(double dt, int test_int) {
//    if (test_int % 2 == 1) {
//        for (auto & planet : planets) {
//
//        planet.x_pos += planet.x_vel * dt;
//        planet.y_pos += planet.y_vel * dt;
//        planet.z_pos += planet.z_vel * dt;
//            //bereken nieuwe versnellingen
//
//
//
//        }
//        }
//    else{
//
//        //bereken de nieuwe versnellingen
//        for (int l = 0; l< planets.size(); l++) {
//            auto &planet = planets[l];
//            std::vector<double> force = {0, 0, 0};
//            //For loop berekent de totale kracht,variabele force[i], die planeet i voelt als gevolg van alle andere planeten
//            for (int j = 0; j < planets.size(); j++) {
//                if (l != j) {
//                    long double dx = planets[j].x_pos - planet.x_pos;
//                    long double dy = planets[j].y_pos - planet.y_pos;
//                    long double dz = planets[j].z_pos - planet.z_pos;
//                    long double r = sqrt(dx * dx + dy * dy + dz * dz);
//                    long double f = G * planet.mass * planets[j].mass / (r * r);
//                    force[0] += f * dx / r;
//                    force[1] += f * dy / r;
//                    force[2] += f * dz / r;
//                    std::cout<<"afstand debug leapfrog "<<r<<std::endl;
//                    std::cout<<"kracht debug "<<force[0]<<std::endl;
//
//                }
//                std::cout<<"acc working debug "<<force[0]/planet.mass<<std::endl;
//
//
//
//            }
//
//            planet.x_acc= force[0] / planet.mass;
//            planet.y_acc= force[1] / planet.mass;
//            planet.z_acc= force[2] / planet.mass;
//
//        }
//        for (auto &planet: planets) {
//
//            planet.x_pos += planet.x_vel * dt;
//            planet.y_pos += planet.y_vel * dt;
//            planet.z_pos += planet.z_vel * dt;
//        }
//
//
//    }
//}

//LEAPFROG TEST
//de richardson test fout is veel beter maar de versnelling klopt niet helemaal

//19 dec 21:36 integraotr

//void Universe::step(double dt, int test_int){
////update positie eerst
//    for (auto & planet : planets) {
//
//        planet.x_pos += planet.x_vel * dt;
//        planet.y_pos += planet.y_vel * dt;
//        planet.z_pos += planet.z_vel * dt;
//    }
//    //bereken nieuwe versnellingen
//    for (int planet_index= 0; planet_index < planets.size(); planet_index++) {
//        auto &planet = planets[planet_index];
//        std::vector<double> force = {0.0, 0.0, 0.0};
//        //For loop berekent de totale kracht,variabele force[i], die planeet i voelt als gevolg van alle andere planeten
//        for (int planet_index_1 = 0; planet_index_1 < planets.size(); planet_index_1++) {
//            if (planet_index != planet_index_1) {
//                double dx = planets[planet_index_1].x_pos - planet.x_pos;
//                double dy = planets[planet_index_1].y_pos - planet.y_pos;
//                double dz = planets[planet_index_1].z_pos - planet.z_pos;
//                double r = sqrt(dx * dx + dy * dy + dz * dz);
//                double f = G * planet.mass * planets[planet_index_1].mass / (r * r);
//                force[0] += f * dx / r;
//                force[1] += f * dy / r;
//                force[2] += f * dz / r;
////                std::cout<<setprecision(15)<<"afstand debug x leapfrog "<<dx<<std::endl;
////                std::cout<<setprecision(15)<<"afstand debug y leapfrog "<<dy<<std::endl;
////                std::cout<<setprecision(15)<<"afstand debug z leapfrog "<<dz<<std::endl;
////                std::cout<<setprecision(15)<<"afstand debug leapfrog "<<r<<std::endl;
////                std::cout<<setprecision(15)<<"kracht debug "<<f<<std::endl;
//            }
//        }
//
//        //de fout zit in dat force[0] niet overeenkomt met de kracht f die planeet 0 voelt
////        std::cout<<"kracht"<<force[0]<<std::endl;
//
//        //leapfrog\
//        leapfrog is niet goed geimplementeerd want je pakt de versnelling niet op het juiste moment
//
//        //vanuit gaande dat de snelheden van planet steeds op een 0.5 te vinden zijn
//        planet.x_acc= force[0] / planet.mass;
//        planet.y_acc= force[1] / planet.mass;
//        planet.z_acc= force[2] / planet.mass;
//
//        //update de snelheden met de nieuwe versnellingen
//        planet.x_vel += planet.x_acc * dt;
//        planet.y_vel += planet.y_acc * dt;
//        planet.z_vel += planet.z_acc * dt;
//
//
//
//        // Met berekende kracht snelheid updaten
//
//    }
//
//}


void Universe::step(double dt, int test_int){
//update positie eerst
    if (test_int%2==0) {
//        std::cout<<"positie update"<<std::endl;
        for (auto &planet: planets) {

            planet.x_pos += planet.x_vel * dt;
            planet.y_pos += planet.y_vel * dt;
            planet.z_pos += planet.z_vel * dt;
        }
        //bereken nieuwe versnellingen
        for (int planet_index = 0; planet_index < planets.size(); planet_index++) {
            auto &planet = planets[planet_index];
            std::vector<double> force = {0.0, 0.0, 0.0};
            //For loop berekent de totale kracht,variabele force[i], die planeet i voelt als gevolg van alle andere planeten
            for (int planet_index_1 = 0; planet_index_1 < planets.size(); planet_index_1++) {
                if (planet_index != planet_index_1) {
                    double dx = planets[planet_index_1].x_pos - planet.x_pos;
                    double dy = planets[planet_index_1].y_pos - planet.y_pos;
                    double dz = planets[planet_index_1].z_pos - planet.z_pos;
                    double r = sqrt(dx * dx + dy * dy + dz * dz);
                    double f = G * planet.mass * planets[planet_index_1].mass / (r * r);
//                    if (planet_index!=2 and planet_index_1!=1){
//                    force[0] += f * dx / r;
//                    force[1] += f * dy / r;
//                    force[2] += f * dz / r;}

                    force[0] += f * dx / r;
                    force[1] += f * dy / r;
                    force[2] += f * dz / r;
//                std::cout<<setprecision(15)<<"afstand debug x leapfrog "<<dx<<std::endl;
//                std::cout<<setprecision(15)<<"afstand debug y leapfrog "<<dy<<std::endl;
//                std::cout<<setprecision(15)<<"afstand debug z leapfrog "<<dz<<std::endl;
//                std::cout<<setprecision(15)<<"afstand debug leapfrog "<<r<<std::endl;
//                std::cout<<setprecision(15)<<"kracht debug "<<f<<std::endl;
                }
            }

            //de fout zit in dat force[0] niet overeenkomt met de kracht f die planeet 0 voelt
//        std::cout<<"kracht"<<force[0]<<std::endl;

            //leapfrog\
        leapfrog is niet goed geimplementeerd want je pakt de versnelling niet op het juiste moment

            //vanuit gaande dat de snelheden van planet steeds op een 0.5 te vinden zijn
            planet.x_acc = force[0] / planet.mass;
            planet.y_acc = force[1] / planet.mass;
            planet.z_acc = force[2] / planet.mass;

        }
    }
    else {
//        std::cout<<"snelheids update"<<std::endl;
        for (auto &planet: planets) {

            planet.x_vel += planet.x_acc * dt;
            planet.y_vel += planet.y_acc * dt;
            planet.z_vel += planet.z_acc * dt;
        }
            //update de snelheden met de nieuwe versnellingen

    }

}


//acceleratie berekenen functie
void Universe::acc_calculate(){
//update positie eerst

//        std::cout<<"positie update"<<std::endl;

        //bereken nieuwe versnellingen
        for (int planet_index = 0; planet_index < planets.size(); planet_index++) {
            auto &planet = planets[planet_index];
            std::vector<double> force = {0.0, 0.0, 0.0};
            //For loop berekent de totale kracht,variabele force[i], die planeet i voelt als gevolg van alle andere planeten
            for (int planet_index_1 = 0; planet_index_1 < planets.size(); planet_index_1++) {
                if (planet_index != planet_index_1) {
                    double dx = planets[planet_index_1].x_pos - planet.x_pos;
                    double dy = planets[planet_index_1].y_pos - planet.y_pos;
                    double dz = planets[planet_index_1].z_pos - planet.z_pos;
                    double r = sqrt(dx * dx + dy * dy + dz * dz);
                    double f = G * planet.mass * planets[planet_index_1].mass / (r * r);


                    force[0] += f * dx / r;
                    force[1] += f * dy / r;
                    force[2] += f * dz / r;
//                std::cout<<setprecision(15)<<"afstand debug x leapfrog "<<dx<<std::endl;
//                std::cout<<setprecision(15)<<"afstand debug y leapfrog "<<dy<<std::endl;
//                std::cout<<setprecision(15)<<"afstand debug z leapfrog "<<dz<<std::endl;
//                std::cout<<setprecision(15)<<"afstand debug leapfrog "<<r<<std::endl;
//                std::cout<<setprecision(15)<<"kracht debug "<<f<<std::endl;
                }
            }
            //vanuit gaande dat de snelheden van planet steeds op een 0.5 te vinden zijn
            planet.x_acc = force[0] / planet.mass;
            planet.y_acc = force[1] / planet.mass;
            planet.z_acc = force[2] / planet.mass;
        }
    }

    //acceleratie berekenen functie

void Universe::vel_calculate(double dt){
        //bereken nieuwe versnellingen
        for (int planet_index = 0; planet_index < planets.size(); planet_index++){
            auto &planet = planets[planet_index];
            //vanuit gaande dat de snelheden van planet steeds op een 0.5 te vinden zijn
            planet.x_vel +=planet.x_acc*dt*0.5;
            planet.y_vel +=planet.y_acc*dt*0.5;
            planet.z_vel +=planet.z_acc*dt*0.5;
        }
    }

void Universe::finalize_velocity(double dt){
    //bereken nieuwe versnellingen
    acc_calculate();
    for (int planet_index = 0; planet_index < planets.size(); planet_index++){
        auto &planet = planets[planet_index];
        //vanuit gaande dat de snelheden van planet steeds op een 0.5 te vinden zijn
        planet.x_vel +=planet.x_acc*dt*0.5;
        planet.y_vel +=planet.y_acc*dt*0.5;
        planet.z_vel +=planet.z_acc*dt*0.5;
    }
}






//
//// define the 4th-order integration constants
//double alpha = 1/(4-pow(2, 4/3));
//double beta_bruh= 1/(2-pow(2, 1/3));
//double gamma_bruh = (1-pow(2, 1/3))/(4-pow(2, 4/3));
//double delta = pow(2, 1/3)/(2-pow(2, 1/3));
//// put them in a vector for "for-loop"
//vector<double> pos_const = {alpha, gamma_bruh, gamma_bruh, alpha};
//vector<double> vel_const = {beta_bruh, delta, beta_bruh};
//
//void Universe::step(double dt, int test_int){
////update positie eerst
//
//    for (int k = 0; k<3; k++) {
//        for (auto &planet: planets) {
//
//            planet.x_pos += pos_const[k] * planet.x_vel * dt;
//            planet.y_pos += pos_const[k] * planet.y_vel * dt;
//            planet.z_pos += pos_const[k] * planet.z_vel * dt;
//        }
//        //bereken nieuwe versnellingen
//        for (int i = 0; i < planets.size(); i++) {
//            auto &planet = planets[i];
//            std::vector<double> force = {0.0, 0.0, 0.0};
//            //For loop berekent de totale kracht,variabele force[i], die planeet i voelt als gevolg van alle andere planeten
//            for (int j = 0; j < planets.size(); j++) {
//                if (i != j) {
//                    double dx = planets[j].x_pos - planet.x_pos;
//                    double dy = planets[j].y_pos - planet.y_pos;
//                    double dz = planets[j].z_pos - planet.z_pos;
//                    double r = sqrt(dx * dx + dy * dy + dz * dz);
//                    double f = G * planet.mass * planets[j].mass / (r * r);
//                    force[0] += f * dx / r;
//                    force[1] += f * dy / r;
//                    force[2] += f * dz / r;
////                    std::cout << setprecision(15) << "afstand debug x leapfrog " << dx << std::endl;
////                    std::cout << setprecision(15) << "afstand debug y leapfrog " << dy << std::endl;
////                    std::cout << setprecision(15) << "afstand debug z leapfrog " << dz << std::endl;
////                    std::cout << setprecision(15) << "afstand debug leapfrog " << r << std::endl;
////                    std::cout << setprecision(15) << "kracht debug " << f << std::endl;
//                }
//            }
//
//            //de fout zit in dat force[0] niet overeenkomt met de kracht f die planeet 0 voelt
//            std::cout << "kracht" << force[0] << std::endl;
//
//            //leapfrog\
//            leapfrog is niet goed geimplementeerd want je pakt de versnelling niet op het juiste moment
//
//            //vanuit gaande dat de snelheden van planet steeds op een 0.5 te vinden zijn
//            planet.x_acc = force[0] / planet.mass;
//            planet.y_acc = force[1] / planet.mass;
//            planet.z_acc = force[2] / planet.mass;
//
//            //update de snelheden met de nieuwe versnellingen
//            planet.x_vel += vel_const[k] * planet.x_acc * dt;
//            planet.y_vel += vel_const[k] * planet.y_acc * dt;
//            planet.z_vel += vel_const[k] * planet.z_acc * dt;
//
//
//
//            // Met berekende kracht snelheid updaten
//
//        }
//    }
//    for (auto &planet: planets) {
//        planet.x_pos += pos_const[3] * planet.x_vel * dt;
//        planet.y_pos += pos_const[3] * planet.y_vel * dt;
//        planet.z_pos += pos_const[3] * planet.z_vel * dt;
//    }
//    }

////Forward euler correct
//void Universe::step(double dt, int test_int) {
//    for (int i = 0; i < planets.size(); i++) {
//        auto &planet = planets[i];
//        std::vector<double> force = {0.0, 0.0, 0.0};
//        //For loop berekent de totale kracht,variabele force[i], die planeet i voelt als gevolg van alle andere planeten
//        for (int j = 0; j < planets.size(); j++) {
//            if (i != j) {
//                double dx = planets[j].x_pos - planet.x_pos;
//                double dy = planets[j].y_pos - planet.y_pos;
//                double dz = planets[j].z_pos - planet.z_pos;
//                double r = sqrt(dx * dx + dy * dy + dz * dz);
//                double f = G * planet.mass * planets[j].mass / (r * r);
//                force[0] += f * dx / r;
//                force[1] += f * dy / r;
//                force[2] += f * dz / r;
//            }
//        }
//
//        //leapfrog\
//        leapfrog is niet goed geimplementeerd want je pakt de versnelling niet op het juiste moment
//
//        //vanuit gaande dat de snelheden van planet steeds op een 0.5 te vinden zijn
//        planet.x_acc= force[0] / planet.mass;
//        planet.y_acc= force[1] / planet.mass;
//        planet.z_acc= force[2] / planet.mass;
//
//
//        planet.x_vel += planet.x_acc * dt;
//        planet.y_vel += planet.y_acc * dt;
//        planet.z_vel += planet.z_acc * dt;
//
//
//
//        // Met berekende kracht snelheid updaten
//
//    }
//    for (auto & planet : planets) {
//
//        planet.x_pos += planet.x_vel * dt;
//        planet.y_pos += planet.y_vel * dt;
//        planet.z_pos += planet.z_vel * dt;
//            //bereken nieuwe versnellingen
//
//
//
//        }
//    }

//void Universe::step(double dt, int test_int) {
//    for (int i = 0; i < planets.size(); i++) {
//        auto &planet = planets[i];
//        std::vector<double> force = {0, 0, 0};
//        //For loop berekent de totale kracht,variabele force[i], die planeet i voelt als gevolg van alle andere planeten
//        for (int j = 0; j < planets.size(); j++) {
//            if (i != j) {
//                double dx = planets[j].x_pos - planet.x_pos;
//                double dy = planets[j].y_pos - planet.y_pos;
//                double dz = planets[j].z_pos - planet.z_pos;
//                double r = sqrt(dx * dx + dy * dy + dz * dz);
//                double f = G * planet.mass * planets[j].mass / (r * r);
//                force[0] += f * dx / r;
//                force[1] += f * dy / r;
//                force[2] += f * dz / r;
//            }
//        }
//
//        //leapfrog\
//        leapfrog is niet goed geimplementeerd want je pakt de versnelling niet op het juiste moment
//
//        //vanuit gaande dat de snelheden van planet steeds op een 0.5 te vinden zijn
//
//
//        planet.x_vel += (force[0] / planet.mass) * dt;
//        planet.y_vel += (force[1] / planet.mass) * dt;
//        planet.z_vel += (force[2] / planet.mass) * dt;
//
//        // Met berekende kracht snelheid updaten
//
//    }
//    for (auto & planet : planets) {
//
//        planet.x_pos += planet.x_vel * dt;
//        planet.y_pos += planet.y_vel * dt;
//        planet.z_pos += planet.z_vel * dt;
//        //bereken nieuwe versnellingen
//
//
//
//    }
//}