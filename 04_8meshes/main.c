# include <stdio.h>
# include <math.h>

int main(void){

    double F = -4.0;

    /*節点の座標を設定*/
    double node_coord_x[9] = {0.0,0.5,1.0,0.0,0.5,1.0,0.0,0.5,1.0};
    double node_coord_y[9] = {0.0,0.0,0.0,0.5,0.5,0.5,1.0,1.0,1.0};

    /*要素と節点の関係を設定*/
    int cnn[8][3] = {
        {0,4,3},{0,1,4},{1,5,4},{1,2,5},{3,7,6},{3,4,7},{4,8,7},{4,5,8}
    };

    /*境界節点の識別*/
    int boarder_node[] = {0,1,2,3,5,6,7,8};
    int non_boarder_node[] = {4};
    int n = 8;
    int m = 1;

    /*解*/
    double u[9];
    for (int i=0; i<9; i++){
        double x = node_coord_x[i];
        double y = node_coord_y[i];
        if (x==0) u[i] = y*y;
        else if (y==0) u[i] = x*x;
        else if (x==1) u[i] = 1+y*y;
        else if (y==1) u[i] = 1+x*x;
        else u[i] = 0;
    }

    /*設定を表示*/
    printf("MESH DATA***************************************************************************\n");
    for (int elem=0; elem<8; elem++){
        int node0 = cnn[elem][0];
        int node1 = cnn[elem][1];
        int node2 = cnn[elem][2];
        printf("%2d   %2d(%lf, %lf)  %2d(%lf, %lf)  %2d(%lf, %lf)\n",elem,
        node0,node_coord_x[node0],node_coord_y[node0],
        node1,node_coord_x[node1],node_coord_y[node1],
        node2,node_coord_x[node2],node_coord_y[node2]);
    }

    /*拡大係数行列の宣言と初期化*/
    double largeA[9][9];
    for (int i=0; i<9; i++) for (int j=0; j<9; j++) largeA[i][j] = 0;
    /*拡大定数項ベクトルの宣言と初期化*/
    double largef[9];
    for (int i=0; i<9; i++) largef[i] = 0;

    /*拡大係数行列と拡大定数項ベクトルの構成*/
    for (int elem=0; elem<8; elem++){
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
    printf("COEF. MATRIX************************************************************************\n");
    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++) printf("%6.3f ",largeA[i][j]);
        printf("\n");
    }

    /*拡大定数項ベクトルの表示*/
    printf("CONST. VECTOR***********************************************************************\n");
    for (int i=0; i<9; i++) printf("%6.3f ",largef[i]);

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

    u[4] = tmpf[0]/tmpA[0][0];

    printf("\nAPPROXIMATES************************************************************************\n");
    for (int i=2; i>=0; i--){
        for (int j=0; j<3; j++) printf("%llf ", u[3*i+j]);
        printf("\n");
    }

    double err = (node_coord_x[4]*node_coord_x[4]+node_coord_y[4]*node_coord_y[4]) - u[4];
    printf("%e", err);
}