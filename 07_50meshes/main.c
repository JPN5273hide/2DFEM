# include <stdio.h>
# include <math.h>

int main(void){

    double F = -4.0;

    /*節点の座標を設定*/
    double node_coord_x[36];
    double node_coord_y[36];
    for (int node=0; node<36; node++){
        node_coord_x[node] = (node%6)/5.0;
        node_coord_y[node] = (node/6)/5.0;
    }

    /*要素と節点の関係を設定*/
    int cnn[50][3] = {
        {0,7,6}, {0,1,7}, {1,8,7}, {1,2,8}, {2,9,8}, {2,3,9}, {3,10,9}, {3,4,10}, {4,11,10}, {4,5,11},
        {6,13,12}, {6,7,13}, {7,14,13}, {7,8,14}, {8,15,14}, {8,9,15}, {9,16,15}, {9,10,16}, {10,17,16}, {10,11,17},
        {12,19,18}, {12,13,19}, {13,20,19}, {13,14,20}, {14,21,20}, {14,15,21}, {15,22,21}, {15,16,22}, {16,23,22}, {16,17,23},
        {18,25,24}, {18,19,25}, {19,26,25}, {19,20,26}, {20,27,26}, {20,21,27}, {21,28,27}, {21,22,28}, {22,29,28}, {22,23,29},
        {24,31,30}, {24,25,31}, {25,32,31}, {25,26,32}, {26,33,32}, {26,27,33}, {27,34,33}, {27,28,34}, {28,35,34}, {28,29,35}
    };

    /*非境界節点*/
    int boarder_node[] = {0,1,2,3,4,5,6,11,12,17,18,23,24,29,30,31,32,33,34,35};
    int non_boarder_node[] = {7,8,9,10,13,14,15,16,19,20,21,22,25,26,27,28};
    int n = 20;
    int m = 16;

    /*解*/
    double u[36];
    for (int i=0; i<36; i++){
        double x = node_coord_x[i];
        double y = node_coord_y[i];
        if (x==0) u[i] = y*y;
        else if (y==0) u[i] = x*x;
        else if (x==1) u[i] = 1+y*y;
        else if (y==1) u[i] = 1+x*x;
        else u[i] = 0;
    }

    /*解析解*/
    double uu[36];
    for (int i=0; i<36; i++){
        double x = node_coord_x[i];
        double y = node_coord_y[i];
        uu[i] = x*x+y*y;
    }

    /*設定を表示*/
    printf("MESH DATA***************************************************************************\n");
    for (int elem=0; elem<36; elem++){
        int node0 = cnn[elem][0];
        int node1 = cnn[elem][1];
        int node2 = cnn[elem][2];
        printf("%2d   %2d(%lf, %lf)  %2d(%lf, %lf)  %2d(%lf, %lf)\n",elem,
        node0,node_coord_x[node0],node_coord_y[node0],
        node1,node_coord_x[node1],node_coord_y[node1],
        node2,node_coord_x[node2],node_coord_y[node2]);
    }

    /*拡大係数行列の宣言と初期化*/
    double largeA[36][36];
    for (int i=0; i<36; i++) for (int j=0; j<36; j++) largeA[i][j] = 0;
    /*拡大定数項ベクトルの宣言と初期化*/
    double largef[36];
    for (int i=0; i<36; i++) largef[i] = 0;

    /*拡大係数行列と拡大定数項ベクトルの構成*/
    for (int elem=0; elem<50; elem++){
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
        printf("ELEM%2d COEF. MATRIX*****************************************************************\n",elem);
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) printf("%5.3f  ",A[i][j]);
            printf("\n");
        }

        /*要素係数行列を拡大係数行列へ組み込み*/
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) largeA[node[i]][node[j]] += A[i][j];
        for (int i=0; i<3; i++) largef[node[i]] += f[i];
    }

    /*拡大要素係数行列の表示*/
    printf("COEF. MATRIX************************************************************************\n");
    for (int i=0; i<36; i++){
        for (int j=0; j<36; j++) printf("%6.3f ",largeA[i][j]);
        printf("\n");
    }

    /*拡大定数項ベクトルの表示*/
    printf("CONST. VECTOR***********************************************************************\n");
    for (int i=0; i<36; i++) printf("%6.3f ",largef[i]);

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
    for (int i=5; i>=0; i--){
        for (int j=0; j<6; j++) printf("%lf ", u[6*i+j]);
        printf("\n");
    }

    printf("\nSOLUTIONS***************************************************************************\n");
    for (int i=5; i>=0; i--){
        for (int j=0; j<6; j++) printf("%lf ", uu[6*i+j]);
        printf("\n");
    }

    double err[36];
    for (int i=0; i<36; i++){
        err[i] = uu[i] - u[i];
    }

    printf("APPROXIMATION ERROR*****************************************************************\n");
    for (int i=5; i>=0; i--){
        for (int j=0; j<6; j++) printf("%e ", err[6*i+j]);
        printf("\n");
    }

    double e = 0.5 - ((u[14]+u[21])/2);
    double f = (u[14]+u[21])/2;
    printf("%e %f", e, f);
}