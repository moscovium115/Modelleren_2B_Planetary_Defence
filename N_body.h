    //
    // Created by anis on 30-11-22.
    //
    #include <vector>
    #ifndef N_BODY_N_BODY_H
    #define N_BODY_N_BODY_H

    class Planet{
    public:
        double mass;
        double radius;
        double x_pos;
        double y_pos;
        double z_pos;
        double x_vel;
        double y_vel;
        double z_vel;
        double x_acc;
        double y_acc;
        double z_acc;
        double x_vel_old;
        double y_vel_old;
        double z_vel_old;
    };

    class Universe {
    public:
        //Attribute planets vector van planeten
        std::vector<Planet> planets;
        //Methode die de tijdintegratie uitvoert
        void step(double dt, int test_int);
        void use_initial_conditions_func(double &t);
        void energies(std::vector<double> &total_impulse,std::vector<double> &ang_momentum, double &sum_kin_energy);
        void detect_collision();
        void detect_collision_2(int &botsing_aantal);
        void inelastic_collision(int i);
        void save_state(double &t);
        void update_dt(double dt, bool dt_changed);
        void acc_calculate();
        void vel_calculate(double dt);
        void finalize_velocity(double dt);
    private:
        long double G = 6.67408e-11;

    };


    #endif //N_BODY_N_BODY_H
