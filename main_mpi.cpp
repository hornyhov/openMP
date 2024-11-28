#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "/opt/homebrew/Cellar/open-mpi/5.0.6/include/mpi.h"
#include <time.h>

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

int main(int argc, char **argv){
    int M = 41, N = 41;
    double A1 = -3.0, A2 = 0.0, B1 = 3.0, B2 = 3.0;
    double EPS;
    double ACC = 1e-6;
    double h1 = (B1 - A1) / M;
    double h2 = (B2 - A2) / N;
    int i, j;
    if (h1 > h2) {
        EPS = h1 * h1;
    }
    else {
        EPS = h2 * h2;
    }

    int rank, size;

    clock_t start_time = clock();

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int k = size / 2 + size % 2;
    int n = size / k;
    int obj_i = rank / k;
    int obj_j = rank % k;

    int row_start = N / n * obj_i;
    int row_end = N / n + 1 + obj_i * (N / n + 1) - n * obj_i;
    int col_start = N / k * obj_j;
    int col_end = N / k + 1 + obj_j * (N / k + 1) - k * obj_j;

    if(row_end > N){
        row_end = N - 1;
    }
    if(col_end > N){
        col_end = N -  1;
    }

    double *data1 = (double*) malloc(M * N * sizeof(double));
    double **F = (double**) malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++)
        F[i] = & (data1[i * N]);
    double *data2 = (double*) malloc (M * N * sizeof(double));
    double **a = (double**) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        a[i] = & (data2[i * N]);    
    double *data3 = (double*) malloc (M * N * sizeof(double));
    double **b = (double**) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        b[i] = & (data3[i * N]);
    double *data4 = (double*) malloc (M * N * sizeof(double));
    double **Ar = (double**) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        Ar[i] = & (data4[i * N]);
    double *data5 = (double*) malloc (M * N * sizeof(double));
    double **r = (double**) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        r[i] = & (data5[i * N]);
    double *data6 = (double*) malloc (M * N * sizeof(double));
    double **w = (double**) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        w[i] = & (data6[i * N]);

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            w[i][j] = 0.0;
            r[i][j] = 0.0;
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            F[i][j] = 0.0;
            Ar[i][j] = 0.0;
        }
    }
    
    double xi, yj, lij, gij, tmp;
    for (i = 1; i < M; i++) {
        for (j = 1; j < N; j++) {
            if(i == 0 || j == 0){
                continue;
            }
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
    double global_norm, global_norm_dr, global_norm_r;
    int iter = 0, msg_iter = 0;

    double *right_send_buf, *right_recv_buf, *left_send_buf, *left_recv_buf, *top_send_buf, *top_recv_buf, *down_send_buf, *down_recv_buf;

    if(obj_j + 1 < k){
        right_send_buf = (double *) malloc((row_end - row_start - 1) * sizeof(double));
        right_recv_buf = (double *) malloc((row_end - row_start - 1) * sizeof(double));
    }
    if(obj_j - 1 >= 0){
        left_send_buf = (double *) malloc((row_end - row_start - 1) * sizeof(double));
        left_recv_buf = (double *) malloc((row_end - row_start - 1) * sizeof(double));
    }
    if(obj_i - 1 >= 0){
        top_send_buf = (double *) malloc((col_end - col_start - 1) * sizeof(double));
        top_recv_buf = (double *) malloc((col_end - col_start - 1) * sizeof(double));
    }
    if(obj_i + 1 < n){
        down_send_buf = (double *) malloc((col_end - col_start - 1) * sizeof(double));
        down_recv_buf = (double *) malloc((col_end - col_start - 1) * sizeof(double));
    }

    bool init_flag = true;
    while (global_norm > ACC || init_flag) {
        init_flag = false;
        norm_r = 0.0, norm_dr = 0.0, norm = 0.0, global_norm = 0.0, global_norm_dr = 0.0, global_norm_r = 0.0;
        MPI_Request recv_right_r, recv_left_r, recv_top_r, recv_down_r, recv_right_Ar, recv_left_Ar, recv_top_Ar, recv_down_Ar;
        int right_send_iter = 0, left_send_iter = 0, top_send_iter = 0, down_send_iter = 0, right_recv_iter = 0, left_recv_iter = 0, top_recv_iter = 0, down_recv_iter = 0;
        if(obj_j + 1 < k){
            MPI_Irecv(right_recv_buf, row_end - row_start - 1, MPI_DOUBLE, obj_i * k + (obj_j + 1) % n, msg_iter, MPI_COMM_WORLD, &recv_right_r);
        }
        if(obj_j - 1 >= 0){
            MPI_Irecv(left_recv_buf, row_end - row_start - 1, MPI_DOUBLE, obj_i * k + (obj_j - 1) % n, msg_iter, MPI_COMM_WORLD, &recv_left_r);
        }
        if(obj_i - 1 >= 0){
            MPI_Irecv(top_recv_buf, col_end - col_start - 1, MPI_DOUBLE, (obj_i - 1) * k + obj_j % n, msg_iter, MPI_COMM_WORLD, &recv_top_r);
        }
        if(obj_i + 1 < n){
            MPI_Irecv(down_recv_buf, col_end - col_start - 1, MPI_DOUBLE, (obj_i + 1) * k + obj_j % n, msg_iter, MPI_COMM_WORLD, &recv_down_r);
        }

        for (i = row_start + 1; i < row_end; i++) {
            for (j = col_start + 1; j < col_end; j++) {
                if(i == 0 || j == 0 || i == M-1 || j == N-1){
                    continue;
                }
                r[i][j] = -(a[i+1][j] * (w[i+1][j] - w[i][j]) / h1 - a[i][j] * (w[i][j] - w[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (w[i][j+1] - w[i][j]) / h2 - b[i][j] * (w[i][j] - w[i][j-1]) / h2) / h2
                        - F[i][j];
                norm_r += r[i][j] * r[i][j] * h1 * h2;
                if(i == row_start + 1 && obj_i - 1 >= 0){
                    top_send_buf[top_send_iter++] = r[i][j];
                }
                if(i == row_end - 1 && obj_i + 1 < n){
                    down_send_buf[down_send_iter++] = r[i][j];
                }
                if(j == col_start + 1 && obj_j - 1 >= 0){
                    left_send_buf[left_send_iter++] = r[i][j];
                }
                if(j == col_end - 1 && obj_j + 1 < k){
                    right_send_buf[right_send_iter++] = r[i][j];
                }
            }
        }
        if(obj_j - 1 >= 0){
            MPI_Request send_left;
            MPI_Isend(left_send_buf, row_end - row_start - 1, MPI_DOUBLE, obj_i * k + (obj_j - 1) % n, msg_iter, MPI_COMM_WORLD, &send_left);
            MPI_Wait(&recv_left_r, MPI_STATUS_IGNORE);
            for(i = row_start + 1; i < row_end; i++){
                r[i][col_start] = left_recv_buf[left_recv_iter++];
            }
            MPI_Irecv(left_recv_buf, row_end - row_start - 1, MPI_DOUBLE, obj_i * k + (obj_j - 1) % n, msg_iter+1, MPI_COMM_WORLD, &recv_left_Ar);
        }
        if(obj_j + 1 < k){
            MPI_Request send_right;
            MPI_Isend(right_send_buf, row_end - row_start - 1, MPI_DOUBLE, obj_i * k + (obj_j + 1) % n, msg_iter, MPI_COMM_WORLD, &send_right);
            MPI_Wait(&recv_right_r, MPI_STATUS_IGNORE);
            for(i = row_start + 1; i < row_end; i++){
                r[i][col_end] = right_recv_buf[right_recv_iter++];
            }
            MPI_Irecv(right_recv_buf, row_end - row_start - 1, MPI_DOUBLE, obj_i * k + (obj_j + 1) % n, msg_iter+1, MPI_COMM_WORLD, &recv_right_Ar);
        }
        if(obj_i + 1 < n){
            MPI_Request send_down;
            MPI_Isend(down_send_buf, col_end - col_start - 1, MPI_DOUBLE, (obj_i + 1) * k + obj_j % n, msg_iter, MPI_COMM_WORLD, &send_down);
            MPI_Wait(&recv_down_r, MPI_STATUS_IGNORE);
            for(j = col_start + 1; j < col_end; j++){
                r[row_end][j] = down_recv_buf[down_recv_iter++];
            }
            MPI_Irecv(down_recv_buf, col_end - col_start - 1, MPI_DOUBLE, (obj_i + 1) * k + obj_j % n, msg_iter+1, MPI_COMM_WORLD, &recv_down_Ar);
        }
        if(obj_i - 1 >= 0){
            MPI_Request send_up;
            MPI_Isend(top_send_buf, col_end - col_start - 1, MPI_DOUBLE, (obj_i - 1) * k + obj_j % n, msg_iter, MPI_COMM_WORLD, &send_up);
            MPI_Wait(&recv_top_r, MPI_STATUS_IGNORE);
            for(j = col_start + 1; j < col_end; j++){
                r[row_start][j] = top_recv_buf[top_recv_iter++];
            }
            MPI_Irecv(top_recv_buf, col_end - col_start - 1, MPI_DOUBLE, (obj_i - 1) * k + obj_j % n, msg_iter+1, MPI_COMM_WORLD, &recv_top_Ar);
        }
        left_send_iter = 0, right_send_iter = 0, top_send_iter = 0, down_send_iter = 0, right_recv_iter = 0, left_recv_iter = 0, top_recv_iter = 0, down_recv_iter = 0;
        msg_iter++;
        for (i = row_start + 1; i < row_end; i++) {
            for (j = col_start + 1; j < col_end; j++) {
                if(i == 0 || j == 0 || i == M-1 || j == N-1){
                    continue;
                }
                Ar[i][j] = -(a[i+1][j] * (r[i+1][j] - r[i][j]) / h1 - a[i][j] * (r[i][j] - r[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (r[i][j+1] - r[i][j]) / h2 - b[i][j] * (r[i][j] - r[i][j-1]) / h2) / h2;
                norm_dr += Ar[i][j] * r[i][j] * h1 * h2;
            }
        }

        MPI_Allreduce(&norm_r, &global_norm_r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&norm_dr, &global_norm_dr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        tau = global_norm_r / global_norm_dr;

        for (i = row_start + 1; i < row_end; i++) {
            for (j = col_start + 1; j < col_end; j++) {
                if(i == 0 || j == 0){
                    continue;
                }
                tmp = w[i][j];
                w[i][j] -= tau * r[i][j];
                tmp = w[i][j] - tmp;
                norm += tmp * tmp * h1 * h2;

                if(i == row_start + 1 && obj_i - 1 >= 0){
                    top_send_buf[top_send_iter++] = w[i][j];
                }
                if(i == row_end - 1 && obj_i + 1 < n){
                    down_send_buf[down_send_iter++] = w[i][j];
                }
                if(j == col_start + 1 && obj_j - 1 >= 0){
                    left_send_buf[left_send_iter++] = w[i][j];
                }
                if(j == col_end - 1 && obj_j + 1 < k){
                    right_send_buf[right_send_iter++] = w[i][j];
                }
            }
        }
        if(obj_j + 1 < k){
            MPI_Request send_right;
            MPI_Isend(right_send_buf, row_end - row_start - 1, MPI_DOUBLE, obj_i * k + (obj_j + 1) % n, msg_iter, MPI_COMM_WORLD, &send_right);
            MPI_Wait(&recv_right_Ar, MPI_STATUS_IGNORE);
            for(i = row_start + 1; i < row_end; i++){
                w[i][col_end] = right_recv_buf[right_recv_iter++];
            }
        }
        if(obj_j - 1 >= 0){
            MPI_Request send_left;
            MPI_Isend(left_send_buf, row_end - row_start - 1, MPI_DOUBLE, obj_i * k + (obj_j - 1) % n, msg_iter, MPI_COMM_WORLD, &send_left);
            MPI_Wait(&recv_left_Ar, MPI_STATUS_IGNORE);
            for(i = row_start + 1; i < row_end; i++){
                w[i][col_start] = left_recv_buf[left_recv_iter++];
            }
        }
        if(obj_i - 1 >= 0){
            MPI_Request send_up;
            MPI_Isend(top_send_buf, col_end - col_start - 1, MPI_DOUBLE, (obj_i - 1) * k + obj_j % n, msg_iter, MPI_COMM_WORLD, &send_up);
            MPI_Wait(&recv_top_Ar, MPI_STATUS_IGNORE);
            for(j = col_start + 1; j < col_end; j++){
                w[row_start][j] = top_recv_buf[top_recv_iter++];
            }
        }
        if(obj_i + 1 < n){
            MPI_Request send_down;
            MPI_Isend(down_send_buf, col_end - col_start - 1, MPI_DOUBLE, (obj_i + 1) * k + obj_j % n, msg_iter, MPI_COMM_WORLD, &send_down);
            MPI_Wait(&recv_down_Ar, MPI_STATUS_IGNORE);
            for(j = col_start + 1; j < col_end; j++){
                w[row_end][j] = down_recv_buf[down_recv_iter++];
            }
        }
        MPI_Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_norm = sqrt(global_norm);
        iter++; msg_iter++;
    }

    printf("Iterations: %d\n", iter);
    printf("Time spent: %f\n", (double)(clock() - start_time) / CLOCKS_PER_SEC);
    
    if(!rank){
        FILE* f = fopen("res.csv", "wb+");
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                // printf("%f ", w[i][j]);
                fprintf(f, "%f,", w[i][j]);
            }
            // printf("\n");
            fprintf(f, "\n");
        }
        fclose(f);
    }

    MPI_Finalize();

    return 0;
}