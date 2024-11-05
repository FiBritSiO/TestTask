
#include <iostream>
#include <stdexcept>
#include <vector>
#include <functional>
#include <algorithm>

struct knurbs
{
    double* ctp;        /* Контрольные точки кривой { x,y,z }            */
    int     kctp;       /* Число контрольных точек, д.б. >2              */
    double* knot;       /* Массив кнот. Количество=kctp+deg-1. М.б.=NULL */
    double* wgt;        /* Массив весов контрольных точек.               */
    /* Количество = числу опорных точек.             */
    /* Если = NULL, то сплайн не рациональный        */
    int    deg;         /* Степень кривой                                */
};

struct Point3D
{
    double x, y, z;


    Point3D operator*(double multiplier)
    {
        return { x * multiplier, y * multiplier, z * multiplier };
    }
};

extern double Precsn = 0.001;

double N(int i, unsigned p, double u, double* knots, unsigned knot_size)
{
    if (i >= knot_size) {

        throw std::out_of_range("Out of range error occured !!");
    }

    if (p == 0)
    {
        if (i + 1 >= knot_size) // return 0.0;
            throw std::out_of_range("Out of range error occured !!");


        if (knots[i] <= u && u < knots[i + 1])  return 1.0;
        else return 0.0;


    }
    else
    {
        double left_part_sum = 0.0;
        double right_part_sum = 0.0;

        if (i + p >= knot_size) // return 0.0;
            throw std::out_of_range("Out of range error occured !!");

        if (knots[i + p] != knots[i])
        {
            left_part_sum = (u - knots[i]) / (knots[i + p] - knots[i]) * N(i, p - 1, u, knots, knot_size);
        }

        if (i + p + 1 >= knot_size) // return left_part_sum;
            throw std::out_of_range("Out of range error occured !!");

        if (knots[i + p + 1] != knots[i + 1])
        {
            right_part_sum = (knots[i + p + 1] - u) / (knots[i + p + 1] - knots[i + 1]) * N(i + 1, p - 1, u, knots, knot_size);
        }

        return left_part_sum + right_part_sum;
    }

}



Point3D evaluatePoint(knurbs* kr, double u)
{
    Point3D formula_numerator = { 0.0, 0.0, 0.0 };
    double formula_denominator = 0.0;
    double* weights;
    double* knot;
    bool weights_was_null = false;
    bool knot_was_null = false;
    unsigned knot_size = kr->kctp + kr->deg + 1;
    if (kr->wgt == NULL)
    {
        weights_was_null = true;
        weights = new double[kr->kctp];
        for (int i = 0; i < kr->kctp; i++)
        {
            weights[i] = 1.0;

        }
    }
    else
    {
        weights = kr->wgt;
    }

    if (kr->knot == NULL)
    {
        knot_was_null = true;
        knot = new double[knot_size];
        for (int i = 0; i < kr->deg + 1; i++)
            knot[i] = 0.0;

        double h = 1.0 / (knot_size - 2 * (kr->deg + 1) + 1);
        int k = kr->deg + 1;
        for (; k < knot_size - (kr->deg + 1); k++)
            knot[k] = knot[k - 1] + h;

        for (; k < knot_size; k++)
            knot[k] = 1.0;

      
    }
    else
    {
        knot = kr->knot;
    }


    for (unsigned i = 0; i < kr->kctp; i++)
    {
        double N_w = N(i, kr->deg, u, knot, knot_size) * weights[i];
        formula_numerator.x = formula_numerator.x + N_w * kr->ctp[3 * i];
        formula_numerator.y = formula_numerator.y + N_w * kr->ctp[3 * i + 1];
        formula_numerator.z = formula_numerator.z + N_w * kr->ctp[3 * i + 2];
        formula_denominator += N_w;
    }

    if (weights_was_null) delete[] weights;
    if (knot_was_null) delete[] knot;
    return formula_numerator * (1.0 / formula_denominator);
}

double distancePoint(knurbs* kr, Point3D tz, double u)
{
    Point3D ev_point = evaluatePoint(kr, u);
    double part_x = pow(tz.x - ev_point.x, 2);
    double part_y = pow(tz.y - ev_point.y, 2);
    double part_z = pow(tz.z - ev_point.z, 2);
    return sqrt(part_x + part_y + part_z);
}


double alphaSearch(std::function<double(double, double)> f, std::vector<double> x, std::vector<double> p, double e = 0.001)
{
    double a = -0.1;
    double b = 0.1;
    double alpha1 = a + (3 - sqrt(5)) / 2 * (b - a);
    double alpha2 = a + (sqrt(5) - 1) / 2 * (b - a);
    double f1 = f(x[0] + alpha1 * p[0], x[1] + alpha1 * p[1]);
    double f2 = f(x[0] + alpha2 * p[0], x[1] + alpha2 * p[1]);
    double tau = (sqrt(5) - 1) / 2;
    double eps_n = (b - a) / 2;
    while (eps_n > e)
    {

        if (f1 <= f2)
        {
            b = alpha2;
            alpha2 = alpha1;
            f2 = f1;
            alpha1 = b - tau * (b - a);
            f1 = f(x[0] + alpha1 * p[0], x[1] + alpha1 * p[1]);
        }
        else
        {
            a = alpha1;
            alpha1 = alpha2;
            f1 = f2;
            alpha2 = a + tau * (b - a);
            f2 = f(x[0] + alpha2 * p[0], x[1] + alpha2 * p[1]);
        }
        eps_n = tau * eps_n;

    }
    return (a + b) / 2;
}



std::vector<double>  CycleCoordLaunch(std::function<double(double, double)> f, std::vector<double> x, std::vector<double> x1_bounds, std::vector<double> x2_bounds, double eps, int max_iters = 100)
{
    std::vector<double> res;
    std::vector<double> _x{ 0.0, 0.0 };
    std::vector<std::vector<double>> e;
    e.push_back(std::vector<double> {1.0, 0.0});
    e.push_back(std::vector<double> {0.0, 1.0});
    int j = 1;

    for (int i = 0; i < max_iters; i++)
    {

        double alpha = alphaSearch(f, x, e[j - 1], eps);
        _x[0] = x[0] + alpha * e[j - 1][0];
        if (!(_x[0] > x1_bounds[0] && _x[0] < x1_bounds[1])) break;
        _x[1] = x[1] + alpha * e[j - 1][1];
        if (!(_x[1] > x2_bounds[0] && _x[1] < x2_bounds[1])) break;

        if (j < 2) { x[0] = _x[0]; x[1] = _x[1]; j++; }
        else
        {
            double cond = sqrt(pow(x[0] - _x[0], 2) + pow(x[1] - _x[1], 2));
            if (cond <= eps)
            {
                x[0] = _x[0];
                x[1] = _x[1];
                break;
            }
            else
            {
                x[0] = _x[0];
                x[1] = _x[1];
                j = 1;
            }
        }
    }
    res.push_back(x[0]);
    res.push_back(x[1]);
    res.push_back(f(x[0], x[1]));
    return res;

}


double findMinDistanceFromPoint(knurbs* kr, double* tz, double* tkr = NULL)
{

    Point3D _tz = { tz[0] , tz[1], tz[2] };
    double a = kr->knot[0];
    int knot_size = kr->deg + kr->kctp + 1;
    double b = kr->knot[knot_size - 1];

    double x1 = a + (3 - sqrt(5)) / 2 * (b - a);
    double x2 = a + (sqrt(5) - 1) / 2 * (b - a);
    double f1 = distancePoint(kr, _tz, x1);
    double f2 = distancePoint(kr, _tz, x2);
    double tau = (sqrt(5) - 1) / 2;
    double eps_n = (b - a) / 2;
    while (eps_n > Precsn)
    {

        if (f1 <= f2)
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - tau * (b - a);
            f1 = distancePoint(kr, _tz, x1);


        }
        else
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + tau * (b - a);
            f2 = distancePoint(kr, _tz, x2);

        }
        eps_n = tau * eps_n;

    }
    if (tkr != NULL)
    {
        Point3D temp = evaluatePoint(kr, (a + b) / 2);
        tkr[0] = temp.x;
        tkr[1] = temp.y;
        tkr[2] = temp.z;

    }

    return distancePoint(kr, _tz, (a + b) / 2);

}


int isCrossingTouching(knurbs* kr1, knurbs* kr2, double* tpk = NULL)
{
    int result = 0;
    int kr1_knot_size = kr1->kctp + kr1->deg + 1;
    int kr2_knot_size = kr2->kctp + kr2->deg + 1;
    double kr1_u_left;
    double kr1_u_right;
    double kr2_u_left;
    double kr2_u_right;
    int kr1_search_index = 0;
    int kr2_search_index = 0;
    std::vector<double> kr1_intervals;
    std::vector<double> kr2_intervals;
    std::vector<double> minDistanceList;
    std::vector<Point3D> point_list;
    auto findMinDistance = [kr1, kr2](double u1, double u2)
    {
        Point3D ev_point1 = evaluatePoint(kr1, u1);
        Point3D ev_point2 = evaluatePoint(kr2, u2);
        double part_x = pow(ev_point1.x - ev_point2.x, 2);
        double part_y = pow(ev_point1.y - ev_point2.y, 2);
        double part_z = pow(ev_point1.z - ev_point2.z, 2);
        return sqrt(part_x + part_y + part_z);
    };

    double h1, h2;

    if (kr1->knot != NULL)
    {
        h1 = (kr1->knot[kr1_knot_size - 1] - kr1->knot[0]) / 10;

        for (int i  = 0 ; i< 11 ;i++)
        {
            kr1_intervals.emplace_back(i*h1);
        }
    }

    else
    {
        h1 = 1.0 / 10;
        for (double x = 0; x < 1; x += h1)
        {
            kr1_intervals.emplace_back(x);
        }
    }

    if (kr2->knot != NULL)
    {
        h2 = (kr2->knot[kr2_knot_size - 1] - kr2->knot[0]) / 10;

        for (int i = 0; i < 11 ; i++)
        {
            kr2_intervals.emplace_back(h2* i);
        }
    }

    else
    {
        h2 = 1.0 / 10;
        for (double x = 0; x < 1; x += h2)
        {
            kr2_intervals.emplace_back(x);
        }
    }
 



    for (int i = 0; i < kr1_intervals.size() - 1; i++)
    {
        double kr1_u_left = kr1_intervals[i];
        double kr1_u_right = kr1_intervals[i + 1];
        double kr1_u0 = (kr1_u_left + kr1_u_right) / 2;
        std::vector<double> kr1_bounds = { kr1_u_left, kr1_u_right };

        for (int j = 0; j < kr2_intervals.size() - 1; j++)
        {
            double kr2_u_left = kr2_intervals[j];
            double kr2_u_right = kr2_intervals[j + 1];
            double kr2_u0 = (kr2_u_left + kr2_u_right) / 2;
            std::vector<double> u0{ kr1_u0,  kr2_u0 };
            std::vector<double> kr2_bounds = { kr2_u_left, kr2_u_right };
            std::vector<double> minDistanse = CycleCoordLaunch(findMinDistance, u0, kr1_bounds, kr2_bounds, 1e-8);
            point_list.push_back(evaluatePoint(kr1, minDistanse[0]));
            minDistanceList.push_back(minDistanse[2]);
        }
    }

    for (int i = 0; i < minDistanceList.size(); i++)
    {
        if (minDistanceList[i] < 1e-6)
        {
            if (tpk != NULL)
            {

                tpk[0] = point_list[i].x;
                tpk[1] = point_list[i].y;
                tpk[2] = point_list[i].z;
            }

            return 1;
        }
    }

    return 0;
}


void task1_test()
{
    knurbs kr1;
    double ctr_points[] = { 10.000000, 30.000000, 0.000000,
              20.000000, 50.000000, 0.000000,
              40.000000, 50.000000, 0.000000,
              65.000000, 0.000000,  0.000000,
              73.106602, 20.606602, 0.000000,
              81.213203, 41.213203, 0.000000,
              81.213203, 71.213203, 0.000000,
              66.213203, 91.213203, 0.000000,
              45.000000, 105.00000, 0.000000,
              30.000000, 85.000000, 0.000000 };

    double knot_point[14] = { 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 7.0, 7.0 }; // см. The Burbs Book с.81

    



    kr1.deg = 3;
    kr1.kctp = 10;
    kr1.ctp = ctr_points;
    kr1.knot = knot_point;
    kr1.wgt = NULL;

    double tz[3] = { 90.0, 95.0 ,0.0 };
    double tkr[3] = { NAN, NAN, NAN };
    double res = findMinDistanceFromPoint(&kr1, tz, tkr);
    std::cout << "Task 1" << std::endl;
    std::cout << "Min distance: " << res << std::endl;
    std::cout <<"Min point: " << tkr[0] << " " << tkr[1] << " " << tkr[2] << std::endl;
}


void task2_test()
{
    knurbs kr1, kr2;

    kr1.deg = 3;
    kr1.kctp = 4;
    double ctr_points1[] = { 35.0, -20.0,  0.0,
     35.0 , 15.0,  0.0,
     95.0,  15.0,  0.0,
     25.0,  40.0,  0.0 };

    double knot_point1[8] =
    {
        0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0
    };

    kr1.ctp = ctr_points1;
    kr1.knot = knot_point1;
    kr1.wgt = NULL;


    kr2.deg = 2;
    kr2.kctp = 3;
    double ctr_points2[] = {
         65.0,      0.0,       0.0,
        -9.641016, 0.0,       0.0,
        55.0,      37.320508, 0.0 };
    double knot_point2[6] = { 0.0, 0.0, 0.0, 2.6179938779914944, 2.6179938779914944,  2.6179938779914944 };
    double wgt2[3] = { 1.0, 0.2588190451025207, 1.0 };
    kr2.ctp = ctr_points2;
    kr2.knot = knot_point2;
    kr2.wgt = wgt2;

    double tpk[] = { NAN, NAN, NAN };
    int res = isCrossingTouching(&kr1, &kr2, tpk);
    std::cout << "Task2" << std::endl;
    std::cout <<"CroosingOrTouching: " << res << std::endl;
    std::cout << "Crossing point: " << tpk[0] << " " << tpk[1] << " " << tpk[2] << std::endl;
}


int main()
{
    task1_test();
    task2_test();





}

