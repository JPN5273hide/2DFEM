# include <stdio.h>
# include <math.h>

# define NDIM1 100

int main(void){

    double F = -4.0;

    /*節点の座標を設定*/
    double node_coord_x[49];
    double node_coord_y[49];
    for (int node=0; node<49; node++){
        node_coord_x[node] = (node%7)/6.0;
        node_coord_y[node] = (node/7)/6.0;
    }

    /*要素と節点の関係を設定*/
    int cnn[72][3] = {
        { 0, 8, 7},{ 0, 1, 8},{ 1, 9, 8},{ 1, 2, 9},{ 2,10, 9},{ 2, 3,10},{ 3,11,10},{ 3, 4,11},{ 4,12,11},{ 4, 5,12},{ 5,13,12},{ 5, 6,13},
        { 7,15,14},{ 7, 8,15},{ 8,16,15},{ 8, 9,16},{ 9,17,16},{ 9,10,17},{10,18,17},{10,11,18},{11,19,18},{11,12,19},{12,20,19},{12,13,20},
        {14,22,21},{14,15,22},{15,23,22},{15,16,23},{16,24,23},{16,17,24},{17,25,24},{17,18,25},{18,26,25},{18,19,26},{19,27,26},{19,20,27},
        {21,29,28},{21,22,29},{22,30,29},{22,23,30},{23,31,30},{23,24,31},{24,32,31},{24,25,32},{25,33,32},{25,26,33},{26,34,33},{26,27,34},
        {28,36,35},{28,29,36},{29,37,36},{29,30,37},{30,38,37},{30,31,38},{31,39,38},{31,32,39},{32,40,39},{32,33,40},{33,41,40},{33,34,41},
        {35,43,42},{35,36,43},{36,44,43},{36,37,44},{37,45,44},{37,38,45},{38,46,45},{38,39,46},{39,47,46},{39,40,47},{40,48,47},{40,41,48}
    };

    /*境界節点の識別*/
    int boarder_node[] = {0,1,2,3,4,5,6,7,13,14,20,21,27,28,34,35,41,42,43,44,45,46,47,48};
    int non_boarder_node[] = {8,9,10,11,12,15,16,17,18,19,22,23,24,25,26,29,30,31,32,33,36,37,38,39,40,43,44,45,46,47};
    int n = 24;
    int m = 25;

    /*解*/
    double u[49];
    for (int i=0; i<49; i++){
        double x = node_coord_x[i];
        double y = node_coord_y[i];
        if (x==0) u[i] = y*y;
        else if (y==0) u[i] = x*x;
        else if (x==1) u[i] = 1+y*y;
        else if (y==1) u[i] = 1+x*x;
        else u[i] = 0;
    }

    /*解析解*/
    double uu[49];
    for (int i=0; i<49; i++){
        double x = node_coord_x[i];
        double y = node_coord_y[i];
        uu[i] = x*x+y*y;
    }

    /*設定を表示*/
    printf("MESH DATA***************************************************************************\n");
    for (int elem=0; elem<72; elem++){
        int node0 = cnn[elem][0];
        int node1 = cnn[elem][1];
        int node2 = cnn[elem][2];
        printf("%2d   %2d(%lf, %lf)  %2d(%lf, %lf)  %2d(%lf, %lf)\n",elem,
        node0,node_coord_x[node0],node_coord_y[node0],
        node1,node_coord_x[node1],node_coord_y[node1],
        node2,node_coord_x[node2],node_coord_y[node2]);
    }

    /*拡大係数行列の宣言と初期化*/
    double largeA[49][49];
    for (int i=0; i<49; i++) for (int j=0; j<49; j++) largeA[i][j] = 0;
    /*拡大定数項ベクトルの宣言と初期化*/
    double largef[49];
    for (int i=0; i<49; i++) largef[i] = 0;

    /*拡大係数行列と拡大定数項ベクトルの構成*/
    for (int elem=0; elem<72; elem++){
        /*座標データの読み込み*/
        int node[3];
        for (int i=0; i<3; i++) node[i] = cnn[elem][i];
        double x0 = node_coord_x[node[0]];
        double x1 = node_coord_x[node[1]];
        double x2 = node_coord_x[node[2]];
        double y0 = node_coord_y[node[0]];
        double y1 = node_coord_y[node[1]];
        double y2 = node_coord_y[node[2]];

        /*要素面積の計算*/
        double D = x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1);
        double S = fabs(D)/2;

        /*要素係数行列の構成*/
        double A[3][3];
        double b[3] = {y1-y2, y2-y0, y0-y1};
        double c[3] = {x2-x1, x0-x2, x1-x0};
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) A[i][j] = (b[j]*b[i]+c[j]*c[i])/(4.0*S);

        /*要素定数項ベクトルの構成*/
        double f[3];
        for (int i=0; i<3; i++) f[i] = F*S/3.0;

        /*要素係数行列の表示*/
        /*
        printf("ELEM%2d COEF. MATRIX*****************************************************************\n",elem);
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) printf("%5.3f  ",A[i][j]);
            printf("\n");
        }*/

        /*要素係数行列を拡大係数行列へ組み込み*/
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) largeA[node[i]][node[j]] += A[i][j];
        for (int i=0; i<3; i++) largef[node[i]] += f[i];
    }

    /*拡大要素係数行列の表示*/
    /*
    printf("COEF. MATRIX************************************************************************\n");
    for (int i=0; i<49; i++){
        for (int j=0; j<49; j++) printf("%6.3f ",largeA[i][j]);
        printf("\n");
    }
    */

    /*拡大定数項ベクトルの表示*/
    /*
    printf("CONST. VECTOR***********************************************************************\n");
    for (int i=0; i<49; i++) printf("%6.3f ",largef[i]);
    */

    /*境界を除去*/  
    double tmpA[m][m];
    double tmpf[m];
    for (int i=0; i<m; i++) for (int j=0; j<m; j++) tmpA[i][j] = largeA[non_boarder_node[i]][non_boarder_node[j]];
    for (int i=0; i<m; i++) tmpf[i] = largef[non_boarder_node[i]];
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            tmpf[j] -= largeA[non_boarder_node[j]][boarder_node[i]]*u[boarder_node[i]];
        }
    }

    /*境界を考慮した拡大要素係数行列の表示*/
    printf("\nCOEF. MATRIX WITHOUT BOUNDARY*******************************************************\n");
    for (int i=0; i<m; i++){
        for (int j=0; j<m; j++) printf("%6.3f ",tmpA[i][j]);
        printf("\n");
    }

    /*境界を考慮した拡大定数項ベクトルの表示*/
    printf("CONST. VECTOR WITHOUT BOUNDARY******************************************************\n");
    for (int i=0; i<m; i++) printf("%6.3f ",tmpf[i]);

    /*ガウス消去法によりAu=fを解く*/
    for (int i=0; i<m-1; i++){
        for (int j=i+1; j<m; j++){
            double aa = tmpA[j][i]/tmpA[i][i];
            tmpf[j] -= aa*tmpf[i];
            for (int k=i+1; k<m; k++) tmpA[j][k] -= aa*tmpA[i][k];
        }
    }
    tmpf[m-1] /= tmpA[m-1][m-1];
    for (int i=m-2; i>=0; i--){
        for (int j=i+1; j<m; j++) tmpf[i] -= tmpA[i][j] * tmpf[j];
        tmpf[i] /= tmpA[i][i];
    }

    for (int i=0; i<m; i++){
        u[non_boarder_node[i]] = tmpf[i];
    }

    /*
    printf("\nCOEF. MATRIX************************************************************************\n");
    for (int i=0; i<m; i++){
        for (int j=0; j<m; j++) printf("%6.3f ",tmpA[i][j]);
        printf("\n");
    }
    */

    printf("\nAPPROXIMATES************************************************************************\n");
    for (int i=6; i>=0; i--){
        for (int j=0; j<7; j++) printf("%llf ", u[7*i+j]);
        printf("\n");
    }

    /*補間関数*/
    double interpolation(x, y, elem) {
        int node[3];
        for (int i=0; i<3; i++) node[i] = cnn[elem][i];
        double x0 = node_coord_x[node[0]];
        double x1 = node_coord_x[node[1]];
        double x2 = node_coord_x[node[2]];
        double y0 = node_coord_y[node[0]];
        double y1 = node_coord_y[node[1]];
        double y2 = node_coord_y[node[2]];

        double u0 = u[node[0]];
        return void;
    }
}