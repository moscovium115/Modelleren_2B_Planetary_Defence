//
// Created by Anis el Hachmioui on 28/11/2022.
//

#include <cmath>
#include <iostream>
#include <vector>
#include "N_body.h"
#include "N_planet.cpp"
#include <fstream>
#include "iomanip"
using namespace std;

extern "C" long long python_perturbations(double* perturbations){
    std::cout<<"cock "<<perturbations[0]<<std::endl;
    std::cout<<"cock "<<perturbations[1]<<std::endl;
    std::cout<<"cock "<<perturbations[2]<<std::endl;
    std::cout<<"cock "<<perturbations[3]<<std::endl;
    std::cout<<"cock "<<perturbations[4]<<std::endl;


    return 0;
}

extern "C" long long perturbation_data_final(double* body_1, double* body_2, int Nt, double dt, bool load_state, bool save_state) {
    //initialize universe
    Universe universe;
    //instantiate position, velocity and acceleration vectors
    //initialize planets
//    Time
    double t = 0;

    //Mass of the Moon
    //hiermee de initial conditions van terugrekenen in de tijd gebruiken
    //instantiate position, velocity and acceleration vectors
    //initialize planets


    //instantiate position, velocity and acceleration vectors
    Planet sun;
    sun.mass = 1.989e30;
    sun.x_pos = 0;
    sun.y_pos = 0;
    sun.z_pos = 0;
    sun.x_vel = 0;
    sun.y_vel = 0;
    sun.z_vel = 0;
    //Aarde toevoegen aan het universum
    universe.planets.push_back(sun);
    //Jupiter toevoegen aan het universum
    double AU=1.4717e11;
    Planet Earth;
    Earth.mass = 5.9722e24;
    Earth.x_pos = 0;
    Earth.y_pos = 1*AU;
    Earth.z_pos = 0;
//    Earth.x_vel = dt/2*-0.00593147,;
    Earth.x_vel=-30000;
    Earth.y_vel = 0;
    Earth.z_vel = 0;
    Earth.radius = 6371000;
    //Aarde toevoegen aan het universum
    universe.planets.push_back(Earth);
    Planet asteroid;
    asteroid.mass = 7.8e10;
    asteroid.x_pos = 6370000;
    asteroid.y_pos = AU;
    asteroid.z_pos = 0;
//Snelheid approximation van t=0.5
//    asteroid.x_vel = -9.8226*dt/2-0.005930675*dt/2;
    asteroid.x_vel = -35e3;
    asteroid.y_vel = 0;
    asteroid.z_vel = 0;
//Maan toevoegen aan het universum
    universe.planets.push_back(asteroid);
    Planet Jupiter;
    Jupiter.mass = 1.898e27;
    Jupiter.x_pos = 7.4076e11;
    Jupiter.y_pos = 0;
    Jupiter.z_pos = 0;
    Jupiter.x_vel = 0;
    Jupiter.y_vel = 13070.0;
    Jupiter.z_vel = 0;
    universe.planets.push_back(Jupiter);
    std::cout<<"aarde x vel initial value:"<<Earth.x_vel<<std::endl;


             //Maan toevoegen aan het universum


    //voor leapfrog moeten we 2x zo veel iteraties nemen

    std::cout << "END DAY:" << dt << std::endl;
    double sum_kin_energy = 0.0;
    std::vector<double> total_impulse;
    std::vector<double> ang_momentum;
    int num_bodies = 2;
    int start_index_load;

    if (load_state) {
        universe.use_initial_conditions_func(t);
        start_index_load=1;
    }
    else{
        start_index_load=0;
        //halve dt snelheids approximation voor leapfrog
        universe.acc_calculate();
        universe.vel_calculate(dt);
    }
    int Nt_leapfrog=2*Nt;
    asteroid = universe.planets[2];
    Earth= universe.planets[1];
    //inital values in de arrays zetten
    body_1[0] = Earth.x_pos;
    body_1[1]= Earth.y_pos;
    body_1[2] = Earth.z_pos;
    body_2[0] = asteroid.x_pos;
    body_2[1] = asteroid.y_pos;
    body_2[2] = asteroid.z_pos;



//    std::cout<<"aarde x acc dbugggggggggggggggg "<<universe.planets[2].x_acc<<std::endl;
    std::cout<<"asteroide x waarde "<< asteroid.x_pos<<std::endl;
    std::cout<<"asteroide y waarde "<< asteroid.y_pos<<std::endl;
    std::cout<<"asteroide z waarde "<< asteroid.z_pos<<std::endl;
    for (int j = 0 ;j < Nt_leapfrog; j++) {
        //j geeft aan welke perturbatie simualitie wordt berekend
        //initialize planets
//        gaat ervan uit dat aarde planeet nummer 2 is en asteroide nummer 4,
        universe.detect_collision();
//Je moet bij een positie update zoals j=Nt de dt veranderen en dan niet integereren en bij d volgende stap

//        Voor het terugrekenen
//        if (j==Nt){
//            std::cout<<"dt changed! "<<j<<std::endl;
//
//            dt=-dt;
//        }
//        if (j!=Nt) {
//            universe.step(dt, j);
//        }
        universe.step(dt,j+start_index_load);




        universe.detect_collision();
//        std::cout<<" C++ dt: "<<dt<<std::endl;

        auto Earth = universe.planets[1];
//        auto jupiter = universe.planets[1];
        auto asteroid = universe.planets[2];
//        auto moon = universe.planets[3];
//        auto asteroid = universe.planets[4];


//        std::cout << "dt!: " << dt << std::endl;
//de x,y,z waardes van de hit worden niet in de array gezet
        if (j%2==0 and j!=Nt_leapfrog){

            body_1[0 + (3) * j/2+3] = Earth.x_pos;
            body_1[1 + (3) * j/2+3] = Earth.y_pos;
            body_1[2 + (3) * j/2+3] = Earth.z_pos;
            body_2[0 + (3) * j/2+3] = asteroid.x_pos;
            body_2[1 + (3) * j/2+3] = asteroid.y_pos;
            body_2[2 + (3) * j/2+3] = asteroid.z_pos;
        }

        if (j==Nt_leapfrog){

//            std::cout<<"aard x pos "<<Earth.x_pos<<std::endl;
//            std::cout<<"aard y pos "<<Earth.y_pos<<std::endl;
//            std::cout<<"aard z pos "<<Earth.z_pos<<std::endl;
//            std::cout<<"aard x vel "<<Earth.x_vel<<std::endl;
//            std::cout<<"aard y vel "<<Earth.y_vel<<std::endl;
//            std::cout<<"aard z vel "<<Earth.z_vel<<std::endl;
//            std::cout<<"asteroide x pos "<<asteroid.x_pos<<std::endl;
//            std::cout<<"asteroide y pos "<<asteroid.y_pos<<std::endl;
//            std::cout<<"asteroide z pos "<<asteroid.z_pos<<std::endl;
//            std::cout<<"asteroide x vel "<<asteroid.x_vel<<std::endl;
//            std::cout<<"asteroide y vel "<<asteroid.y_vel<<std::endl;
//            std::cout<<"asteroide z vel "<<asteroid.z_vel<<std::endl;
//            Earth = universe.planets[1];
//            asteroid = universe.planets[2];
//            std::cout<<"aard x vel"<<universe.planets[2].x_vel<<std::endl;
        }
        if(load_state and j==Nt_leapfrog-1){
            universe.finalize_velocity(dt);
//            Earth = universe.planets[1];
            asteroid = universe.planets[2];
        }

//        std::cout<<"asteroide x waarde "<< asteroid.x_pos<<std::endl;
//        std::cout<<"asteroide y waarde "<< asteroid.y_pos<<std::endl;
//        std::cout<<"asteroide z waarde "<< asteroid.z_pos<<std::endl;
//        std::cout<<"asteroide x vel "<<asteroid.x_vel<<std::endl;
//        std::cout<<"asteroide y vel "<<asteroid.y_vel<<std::endl;
//        std::cout<<"asteroide z vel "<<asteroid.z_vel<<std::endl;
//        std::cout<<"asteroide x acc "<<asteroid.x_acc<<std::endl;
//        std::cout<<"asteroide y acc "<<asteroid.y_acc<<std::endl;
//        std::cout<<"asteroide z acc "<<asteroid.z_acc<<std::endl;

//        std::cout<<setprecision(20)<<"Earth x waarde "<< Earth.x_pos<<std::endl;
//        std::cout<<setprecision(20)<<"Earth y waarde "<< Earth.y_pos<<std::endl;
//        std::cout<<setprecision(20)<<"Earth z waarde "<< Earth.z_pos<<std::endl;
//        std::cout<<setprecision(20)<<"Earth x vel "<< Earth.x_vel<<std::endl;
//        std::cout<<setprecision(20)<<"Earth y vel "<< Earth.y_vel<<std::endl;
//        std::cout<<setprecision(20)<<"Earth z vel "<< Earth.z_vel<<std::endl;
        t += dt;
//        std::cout << "python testje: " << j << std::endl;
    }

    if(save_state){
        universe.save_state(t);
    }
    return 0;
}


extern "C" long long perturbations_simulaties(double* body_1, double* body_2, int Nt, double dt, bool load_state, bool save_state, double* perturbations, int num_perturbations, int indexje, double* prob_arr) {
    int num_of_hits=0;

    for (int pertu_iter = 0; pertu_iter < num_perturbations; pertu_iter++) {


        //initialize universe
        Universe universe;
        //instantiate position, velocity and acceleration vectors
        //initialize planets
//    Time
        double t = 0;

        //Mass of the Moon
        //hiermee de initial conditions van terugrekenen in de tijd gebruiken
        //instantiate position, velocity and acceleration vectors
        //initialize planets


        //instantiate position, velocity and acceleration vectors
        Planet sun;
        sun.mass = 1.989e30;
        sun.x_pos = 0;
        sun.y_pos = 0;
        sun.z_pos = 0;
        sun.x_vel = 0;
        sun.y_vel = 0;
        sun.z_vel = 0;
        //Aarde toevoegen aan het universum
        universe.planets.push_back(sun);
        //Jupiter toevoegen aan het universum
        double AU = 1.4717e11;
        Planet Earth;
        Earth.mass = 5.9722e24;
        Earth.x_pos = 0;
        Earth.y_pos = 1 * AU;
        Earth.z_pos = 0;
//    Earth.x_vel = dt/2*-0.00593147,;
        Earth.x_vel = -30000;
        Earth.y_vel = 0;
        Earth.z_vel = 0;
        Earth.radius = 6371000;
        //Aarde toevoegen aan het universum
        universe.planets.push_back(Earth);
        Planet asteroid;
        asteroid.mass = 7.8e10;
        asteroid.x_pos = 6370000;
        asteroid.y_pos = AU;
        asteroid.z_pos = 0;
//Snelheid approximation van t=0.5
//    asteroid.x_vel = -9.8226*dt/2-0.005930675*dt/2;
        asteroid.x_vel = -35e3;
        asteroid.y_vel = 0;
        asteroid.z_vel = 0;
//Maan toevoegen aan het universum
        universe.planets.push_back(asteroid);
        Planet Jupiter;
        Jupiter.mass = 1.898e27;
        Jupiter.x_pos = 7.4076e11;
        Jupiter.y_pos = 0;
        Jupiter.z_pos = 0;
        Jupiter.x_vel = 0;
        Jupiter.y_vel = 13070.0;
        Jupiter.z_vel = 0;
        universe.planets.push_back(Jupiter);
        std::cout << "aarde x vel initial value:" << Earth.x_vel << std::endl;


        //Maan toevoegen aan het universum


        //voor leapfrog moeten we 2x zo veel iteraties nemen

        std::cout << "END DAY:" << dt << std::endl;
        double sum_kin_energy = 0.0;
        std::vector<double> total_impulse;
        std::vector<double> ang_momentum;
        int num_bodies = 2;
        int start_index_load;

        if (load_state) {
            universe.use_initial_conditions_func(t);
            start_index_load = 1;
            universe.planets[2].x_pos += perturbations[0 + 3 * pertu_iter];
            universe.planets[2].y_pos += perturbations[1 + 3 * pertu_iter];
            universe.planets[2].z_pos += perturbations[2 + 3 * pertu_iter];
        } else {
            start_index_load = 0;
            //halve dt snelheids approximation voor leapfrog
            universe.acc_calculate();
            universe.vel_calculate(dt);
        }
        int Nt_leapfrog = 2 * Nt;
        //zodat hun waardes geupdatet zijn
        asteroid = universe.planets[2];
        Earth = universe.planets[1];
        //inital values in de arrays zetten
        body_1[0] = Earth.x_pos;
        body_1[1] = Earth.y_pos;
        body_1[2] = Earth.z_pos;
        body_2[0] = asteroid.x_pos;
        body_2[1] = asteroid.y_pos;
        body_2[2] = asteroid.z_pos;



//    std::cout<<"aarde x acc dbugggggggggggggggg "<<universe.planets[2].x_acc<<std::endl;
//        std::cout << "asteroide x waarde " << asteroid.x_pos << std::endl;
//        std::cout << "asteroide y waarde " << asteroid.y_pos << std::endl;
//        std::cout << "asteroide z waarde " << asteroid.z_pos << std::endl;
        for (int j = 0; j < Nt_leapfrog; j++) {
            //j geeft aan welke perturbatie simualitie wordt berekend
            //initialize planets
//        gaat ervan uit dat aarde planeet nummer 2 is en asteroide nummer 4,
//            universe.detect_collision(num_of_hits);
//Je moet bij een positie update zoals j=Nt de dt veranderen en dan niet integereren en bij d volgende stap

//        Voor het terugrekenen
//        if (j==Nt){
//            std::cout<<"dt changed! "<<j<<std::endl;
//
//            dt=-dt;
//        }
//        if (j!=Nt) {
//            universe.step(dt, j);
//        }
            universe.step(dt, j + start_index_load);


            universe.detect_collision_2(num_of_hits);

//        std::cout<<" C++ dt: "<<dt<<std::endl;

            auto Earth = universe.planets[1];
//        auto jupiter = universe.planets[1];
            auto asteroid = universe.planets[2];
//        auto moon = universe.planets[3];
//        auto asteroid = universe.planets[4];


//        std::cout << "dt!: " << dt << std::endl;
//de x,y,z waardes van de hit worden niet in de array gezet
            if (j % 2 == 0 and j != Nt_leapfrog) {

                body_1[0 + (3) * j / 2 + 3] = Earth.x_pos;
                body_1[1 + (3) * j / 2 + 3] = Earth.y_pos;
                body_1[2 + (3) * j / 2 + 3] = Earth.z_pos;
                body_2[0 + (3) * j / 2 + 3] = asteroid.x_pos;
                body_2[1 + (3) * j / 2 + 3] = asteroid.y_pos;
                body_2[2 + (3) * j / 2 + 3] = asteroid.z_pos;
            }

            if (j == Nt_leapfrog) {

//            std::cout<<"aard x pos "<<Earth.x_pos<<std::endl;
//            std::cout<<"aard y pos "<<Earth.y_pos<<std::endl;
//            std::cout<<"aard z pos "<<Earth.z_pos<<std::endl;
//            std::cout<<"aard x vel "<<Earth.x_vel<<std::endl;
//            std::cout<<"aard y vel "<<Earth.y_vel<<std::endl;
//            std::cout<<"aard z vel "<<Earth.z_vel<<std::endl;
//            std::cout<<"asteroide x pos "<<asteroid.x_pos<<std::endl;
//            std::cout<<"asteroide y pos "<<asteroid.y_pos<<std::endl;
//            std::cout<<"asteroide z pos "<<asteroid.z_pos<<std::endl;
//            std::cout<<"asteroide x vel "<<asteroid.x_vel<<std::endl;
//            std::cout<<"asteroide y vel "<<asteroid.y_vel<<std::endl;
//            std::cout<<"asteroide z vel "<<asteroid.z_vel<<std::endl;
//            Earth = universe.planets[1];
//            asteroid = universe.planets[2];
//            std::cout<<"aard x vel"<<universe.planets[2].x_vel<<std::endl;
            }
            if (load_state and j == Nt_leapfrog - 1) {
                universe.finalize_velocity(dt);
//            Earth = universe.planets[1];
                asteroid = universe.planets[2];
            }

//        std::cout<<"asteroide x waarde "<< asteroid.x_pos<<std::endl;
//        std::cout<<"asteroide y waarde "<< asteroid.y_pos<<std::endl;
//        std::cout<<"asteroide z waarde "<< asteroid.z_pos<<std::endl;
//        std::cout<<"asteroide x vel "<<asteroid.x_vel<<std::endl;
//        std::cout<<"asteroide y vel "<<asteroid.y_vel<<std::endl;
//        std::cout<<"asteroide z vel "<<asteroid.z_vel<<std::endl;
//        std::cout<<"asteroide x acc "<<asteroid.x_acc<<std::endl;
//        std::cout<<"asteroide y acc "<<asteroid.y_acc<<std::endl;
//        std::cout<<"asteroide z acc "<<asteroid.z_acc<<std::endl;

//        std::cout<<setprecision(20)<<"Earth x waarde "<< Earth.x_pos<<std::endl;
//        std::cout<<setprecision(20)<<"Earth y waarde "<< Earth.y_pos<<std::endl;
//        std::cout<<setprecision(20)<<"Earth z waarde "<< Earth.z_pos<<std::endl;
//        std::cout<<setprecision(20)<<"Earth x vel "<< Earth.x_vel<<std::endl;
//        std::cout<<setprecision(20)<<"Earth y vel "<< Earth.y_vel<<std::endl;
//        std::cout<<setprecision(20)<<"Earth z vel "<< Earth.z_vel<<std::endl;
            t += dt;
//        std::cout << "python testje: " << j << std::endl;
        }

        if (save_state) {
            universe.save_state(t);
        }
    }
    float hit_chance=100.0*num_of_hits/num_perturbations;
    std::cout<<setprecision(8)<<" De botsingskans is : "<<hit_chance<<std::endl;
    prob_arr[indexje]=hit_chance;
        return 0;
    }



extern "C" long long perturbation_data_lyapunov(double* body_1, double* body_2, int Nt, double dt, bool load_state, bool save_state,double* init_diff_vector, double* final_diff_vector) {
    //
    std::vector<double> final_diff_vector_1={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

    //initialize universe
    for (int run=0; run<2; run++){

        Universe universe;
        //instantiate position, velocity and acceleration vectors
        //initialize planets
//    Time
        double t = 0;

        //Mass of the Moon
        //hiermee de initial conditions van terugrekenen in de tijd gebruiken
        //instantiate position, velocity and acceleration vectors
        //initialize planets


        //instantiate position, velocity and acceleration vectors
        Planet sun;
        sun.mass = 1.989e30;
        sun.x_pos = 0;
        sun.y_pos = 0;
        sun.z_pos = 0;
        sun.x_vel = 0;
        sun.y_vel = 0;
        sun.z_vel = 0;
        //Aarde toevoegen aan het universum
        universe.planets.push_back(sun);
        //Jupiter toevoegen aan het universum
        Planet Earth;
        Earth.mass = 5.972e24;
        Earth.x_pos = 1.496e11;
        Earth.y_pos = 0;
        Earth.z_pos = 0;
        Earth.x_vel = dt/2*-0.00593147;
        Earth.y_vel = 29780.0;
        Earth.z_vel = 0;
        Earth.radius = 6371000;
        std::cout<<"aarde x pos initial value:"<<Earth.x_pos<<std::endl;
        //Aarde toevoegen aan het universum
        universe.planets.push_back(Earth);
        Planet asteroid;
        asteroid.mass = 1;
        asteroid.x_pos = 1.496e11-6370000+run*init_diff_vector[0];
        asteroid.y_pos = 0+run*init_diff_vector[1];
        asteroid.z_pos = 0+run*init_diff_vector[2];
//Snelheid approximation van t=0.5
        asteroid.x_vel = 0;
        asteroid.y_vel = 35000;
        asteroid.z_vel = 0;
//Maan toevoegen aan het universum
        universe.planets.push_back(asteroid);
        Planet Jupiter;
        Jupiter.mass = 1.898e27;
        Jupiter.x_pos = 5.20336301e11;
        Jupiter.y_pos = 0;
        Jupiter.z_pos = 0;
        Jupiter.x_vel = 0;
        Jupiter.y_vel = 13070;
        Jupiter.z_vel = 0;
        universe.planets.push_back(Jupiter);
        Planet Saturnus;
        Saturnus.mass = 5.683e26;
        Saturnus.x_pos = 9.53707032e11;
        Saturnus.y_pos = 0;
        Saturnus.z_pos = 0;
        Saturnus.x_vel = 0;
        Saturnus.y_vel = 9690;
        Saturnus.z_vel = 0;
        universe.planets.push_back(Saturnus);


        //Maan toevoegen aan het universum


        //voor leapfrog moeten we 2x zo veel iteraties nemen

        std::cout << "END DAY:" << dt << std::endl;
        double sum_kin_energy = 0.0;
        std::vector<double> total_impulse;
        std::vector<double> ang_momentum;
        int num_bodies = 2;

        if (load_state) {
            universe.use_initial_conditions_func(t);
        }
        int Nt_leapfrog=2*Nt;

        for (int j = 0; j < Nt_leapfrog; j++) {
            //j geeft aan welke perturbatie simualitie wordt berekend
            //initialize planets
//        gaat ervan uit dat aarde planeet nummer 2 is en asteroide nummer 4,
            universe.detect_collision();
            universe.step(dt, j);
//        std::cout<<" C++ dt: "<<dt<<std::endl;

            auto Earth = universe.planets[1];
//        auto jupiter = universe.planets[1];
            auto asteroid = universe.planets[2];
            auto Jupiter = universe.planets[3];
            auto Saturnus = universe.planets[4];
//        auto moon = universe.planets[3];
//        auto asteroid = universe.planets[4];


//        std::cout << "dt!: " << dt << std::endl;
            if (run==1) {
                if (j % 2 == 0) {
                    body_1[0 + (3) * j / 2] = Earth.x_pos;
                    body_1[1 + (3) * j / 2] = Earth.y_pos;
                    body_1[2 + (3) * j / 2] = Earth.z_pos;
                    body_2[0 + (3) * j / 2] = asteroid.x_pos;
                    body_2[1 + (3) * j / 2] = asteroid.y_pos;
                    body_2[2 + (3) * j / 2] = asteroid.z_pos;

                }
            }
            if (j%(Nt_leapfrog/2)){
                dt=dt;
            }

//        std::cout<<"asteroide x waarde "<< asteroid.x_pos<<std::endl;
//        std::cout<<"asteroide y waarde "<< asteroid.y_pos<<std::endl;
//        std::cout<<"asteroide z waarde "<< asteroid.z_pos<<std::endl;
//        std::cout<<"Earth x waarde "<< Earth.x_pos<<std::endl;
//        std::cout<<"Earth y waarde "<< Earth.y_pos<<std::endl;
//        std::cout<<"Earth z waarde "<< Earth.z_pos<<std::endl;









            t += dt;
//        std::cout << "python testje" << j << std::endl;
        }
        final_diff_vector_1[0+run*6]=asteroid.x_pos;
        final_diff_vector_1[1+run*6]=asteroid.y_pos;
        final_diff_vector_1[2+run*6]=asteroid.z_pos;
        final_diff_vector_1[3+run*6]=asteroid.x_vel;
        final_diff_vector_1[4+run*6]=asteroid.y_vel;
        final_diff_vector_1[5+run*6]=asteroid.z_vel;

    }
    std::cout<<"final diff vector 1"<<std::endl;
    print_double(final_diff_vector_1);
    final_diff_vector[0]=final_diff_vector_1[0]-final_diff_vector_1[6];
    final_diff_vector[1]=final_diff_vector_1[1]-final_diff_vector_1[7];
    final_diff_vector[2]=final_diff_vector_1[2]-final_diff_vector_1[8];
    final_diff_vector[3]=final_diff_vector_1[3]-final_diff_vector_1[9];
    final_diff_vector[4]=final_diff_vector_1[4]-final_diff_vector_1[10];
    final_diff_vector[5]=final_diff_vector_1[5]-final_diff_vector_1[11];


    std::cout<<"final diff vector"<<std::endl;
    std::cout<<final_diff_vector[0]<<std::endl;
    std::cout<<final_diff_vector[1]<<std::endl;
    std::cout<<final_diff_vector[2]<<std::endl;


    return 0;
}




