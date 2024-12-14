#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

typedef struct {
    double x, y;
} Point;

typedef struct {
    Point F, G, H, P;
} Rectangle;

typedef struct {
    Point A, B, C, D;
} Trapezoid;

double left_line(Point a){
    return 3.0 * a.x + 9 - a.y;
}

double right_line(Point a){
    return -3.0 * a.x + 9 - a.y;
}

Point calc_left_point(Point a){
    Point res = {(a.y - 9.0) / 3.0, a.y};
    return res;
}

Point calc_right_point(Point a){
    Point res = {(a.y - 9.0) / -3.0, a.y};
    return res;
}

double intersection_area(Rectangle rect) {
    
    Trapezoid trap = {{-3, 0}, {3, 0}, {-2, 3}, {2, 3}};
    
    Point left_point_bottom;
    Point right_point_bottom;

    Point left_point_top;
    Point right_point_top;

    if ((left_line(rect.F) <= 0 && left_line(rect.G) <= 0) || (right_line(rect.F) < 0 && right_line(rect.G) < 0)){
        return 0.0;
    }
    else if (left_line(rect.F) < 0 && left_line(rect.G) > 0 && left_line(rect.H) < 0 && left_line(rect.P) < 0) {
        left_point_top.x = rect.G.x;
        right_point_top.x = rect.G.x;

        left_point_top.y = 3.0 * rect.G.x + 9.0;
        right_point_top.y = 3.0 * rect.G.x + 9.0;

        left_point_bottom = calc_left_point(rect.F);
        right_point_bottom = rect.G;
    }
    else if (right_line(rect.F) > 0 && right_line(rect.G) < 0 && right_line(rect.H) < 0 && right_line(rect.P) < 0) {
        left_point_top.x = rect.F.x;
        right_point_top.x = rect.F.x;

        left_point_top.y = -3.0 * rect.F.x + 9.0;
        right_point_top.y = -3.0 * rect.F.x + 9.0;

        left_point_bottom = rect.F;
        right_point_bottom = calc_right_point(rect.G);
    } else {
        if(left_line(rect.F) >= 0){
            left_point_bottom = rect.F;
        }
        else {
            left_point_bottom = calc_left_point(rect.F);
        }
        if(right_line(rect.G) >= 0){
            right_point_bottom = rect.G;
        }
        else{
            right_point_bottom = calc_right_point(rect.G);
        }
        if(left_line(rect.H) >= 0){
            left_point_top = rect.H;
        }
        else {
            left_point_top = calc_left_point(rect.H);
        }

        if(right_line(rect.P) >= 0){
            right_point_top = rect.P;
        }
        else {
            right_point_top = calc_right_point(rect.P);
        }
    }

    return (right_point_bottom.x - left_point_bottom.x + right_point_top.x - left_point_top.x) / 2.0 * (left_point_top.y - left_point_bottom.y);
}

double lengh_of_vert_line(Point top, Point bottom){

    Point top_point, bottom_point;

    if(left_line(bottom) <= 0 || right_line(bottom) <= 0) {
        return 0.0;
    } else if (left_line(bottom) > 0 && left_line(top) < 0){
        top_point.y = 3.0 * bottom.x + 9.0;
        top_point.x = top.x;
        bottom_point = bottom;
    } else if (right_line(bottom) > 0 && left_line(top) < 0){
        top_point.y = -3.0 * bottom.x + 9.0;
        top_point.x = top.x;
        bottom_point = bottom;
    }
    else {
        top_point = top;
        bottom_point = bottom;
    }

    return (top_point.y - bottom_point.y);
}

double lengh_of_horz_line(Point left, Point right){

    Point left_point, right_point;

    if(left_line(right) <= 0 || right_line(left) <= 0) {
        return 0.0;
    } else if (left_line(right) > 0 && left_line(left) < 0){
        left_point = calc_left_point(right);
        right_point = right;
    } else if (right_line(left) > 0 && left_line(right) < 0){
        left_point = left;
        right_point = calc_right_point(left);
    }
    else {
        left_point = left;
        right_point = right;
    }

    return (right_point.x - left_point.x);
}

double f(double x_l, double y_l, double x_r, double y_r, double h1, double h2) {
    Rectangle rect;

    rect.F = (Point){x_l, y_l};
    rect.G = (Point){x_r, y_l};

    rect.H = (Point){x_l, y_r};
    rect.P = (Point){x_r, y_r};

    return intersection_area(rect) / (h1 * h2);

}

int main() {

    int M = 41,  N = 41;
    double A1 = -3.0, A2 = 0.0, B1 = 3.0, B2 = 3.0;
    double ACC = 1e-6;
    double start_time,  end_time,  EPS;
    double h1 = (B1 - A1) / M;
    double h2 = (B2 - A2) / N;
    if (h1 > h2) {
        EPS = h1 * h1;
    }
    else {
        EPS = h2 * h2;
    }
    double w[M][N], r[M][N], a[M][N], b[M][N], F[M][N], Diff[M][N], Ar[M][N];
    int i, j, k;
    double xi, yj, lij, gij, tmp;

    // start_time = omp_get_wtime();

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            b[i][j] = 0.0; F[i][j] = 0.0; Diff[i][j] = 0.0; Ar[i][j] = 0.0;
            w[i][j] = 0.0; r[i][j] = 0.0; a[i][j] = 0.0;
        }
    }

    for (i = 1; i < M; i++) {
        for (j = 1; j < N; j++) {
            xi = A1 + i*h1;
            yj = A2 + j*h2;

            lij = lengh_of_vert_line((Point){xi-0.5*h1, yj+0.5*h2}, (Point){xi-0.5*h1, yj-0.5*h2});
            gij = lengh_of_horz_line((Point){xi-0.5*h1, yj-0.5*h2}, (Point){xi+0.5*h1, yj-0.5*h2});

            a[i][j] = lij / h2 + (1.0 - lij / h2) / EPS;
            b[i][j] = gij / h1 + (1.0 - gij / h1) / EPS;

            if(i != M-1 && j != N-1) {
                F[i][j] = f(xi-0.5*h1, yj-0.5*h2, xi+0.5*h1, yj+0.5*h2, h1, h2);
            }
        }
    }

    double tau, norm_r, norm_dr, norm;
    int iter = 0;
    bool init_flag = true;

    while (norm > ACC || init_flag) {
        init_flag = false;
        norm_r = 0.0; norm_dr = 0.0; norm = 0.0;

        for (i = 1; i < M-1; i++) {
            for (j = 1; j < N-1; j++) {
                r[i][j] = -(a[i+1][j] * (w[i+1][j] - w[i][j]) / h1 - a[i][j] * (w[i][j] - w[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (w[i][j+1] - w[i][j]) / h2 - b[i][j] * (w[i][j] - w[i][j-1]) / h2) / h2
                        - F[i][j];
                norm_r += r[i][j] * r[i][j] * h1 * h2;
            }
        }

        for (i = 1; i < M-1; i++) {
            for (j = 1; j < N-1; j++) {
                Ar[i][j] = -(a[i+1][j] * (r[i+1][j] - r[i][j]) / h1 - a[i][j] * (r[i][j] - r[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (r[i][j+1] - r[i][j]) / h2 - b[i][j] * (r[i][j] - r[i][j-1]) / h2) / h2;
                norm_dr += Ar[i][j] * r[i][j] * h1 * h2;
            }
        }

        for (i = 1; i < M; i++) {
            for (j = 1; j < N; j++) {
                tmp = w[i][j];
                w[i][j] -= norm_r / norm_dr * r[i][j];
                Diff[i][j] = w[i][j] - tmp;
                norm += Diff[i][j] * Diff[i][j] * h1 * h2;
            }
        }

        ++iter;
        norm = sqrt(norm);
    } 

    // printf("Time spent: %f\n", omp_get_wtime() - start_time);
    printf("Total interations: %d\n", iter);

    std::ofstream myfile;
    myfile.open("res.csv");
    

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            // printf("%f ", w[i][j]);
            myfile << w[i][j] << ",";
        }
        // printf("\n");
        myfile << "\n";
    }
    // myfile.close();

    return 0;
}