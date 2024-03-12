# include <stdio.h>
# include <math.h>

int main(void){

    double F = 4; /*XXX: 多分違う*/
    
    int elem_type[8] = {2,1,2,1,2,1,2,1};
    int cnn[8][3] = {{0,4,3},{0,1,4},{1,5,4},{1,2,5},{3,7,6},{3,4,7},{4,8,7},{4,5,8}}; /*elementとnodeの関係*/

    double node_coord_x[9] = {0.0,0.4,1.0,0.0,0.4,1.0,0.0,0.4,1.0}; /*nodeのx座標*/
    double node_coord_y[9] = {0.0,0.0,0.0,0.7,0.7,0.7,1.0,1.0,1.0}; /*nodeのy座標*/

    double D[8]; /*固有値*/
    double S[8]; /*要素面積*/

    /*HACK: ループ内で定義すればOK*/
    double A[8][3][3]; /*要素係数マトリクス*/
    double f[8][3]; /*要素自由項ベクトル*/

    double largeA[9][9]; /*拡大係数マトリクス*/
    double largef[9]; /*拡大自由項ベクトル*/

    double u[9]; /*解*/

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            largeA[i][j] = 0;
        }
        largef[i] = 0;
    }

    /*各要素に対してA,fを計算 -> 全体のA,fに組み込み*/
    for (int elem=0; elem<8; elem++){
        /*TODO: 直接定義でOK*/
        int node[3];
        for (int i=0; i<3; i++){
            node[i] = cnn[elem][i];
        }

        double x1 = node_coord_x[node[0]];
        double x2 = node_coord_x[node[1]];
        double x3 = node_coord_x[node[2]];
        double y1 = node_coord_y[node[0]];
        double y2 = node_coord_y[node[1]];
        double y3 = node_coord_y[node[2]];

        double b[3] = {y2-y3,y3-y1,y1-y2};
        double c[3] = {x3-x2,x1-x3,x2-x1};

        D[elem] = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
        S[elem] = fabs(D[elem])/2;
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                A[elem][i][j] = (b[j]*b[i]+c[j]*c[i])/(4*S[elem]);
            }
        }
        for (int i=0; i<3; i++){
            f[elem][i] = F*S[elem]/3;
        }

        printf("COEFFICIENT MATRIX for element%d:\n", elem+1);
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                printf("%lf, ",A[elem][i][j]);
            }
            printf("\n");
        }
        printf("CONSTANT TERM VECTOR for element%d:\n%lf\n", elem+1, f[elem][0]);

        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                largeA[node[i]][node[j]] += A[elem][i][j];
            }
            largef[node[i]] += f[elem][i];
        }
    }

    printf("\nCOEFFICIENT MATRIX:\n");
        for (int i=0; i<9; i++){
            for (int j=0; j<9; j++){
                printf("%lf, ",largeA[i][j]);
            }
            printf("\n");
        }
    printf("CONSTANT TERM VECTOR:\n");
    for (int i=0; i<9; i++){
        printf("%lf, ",largef[i]);
    }

    /*基本境界条件の適用*/
    u[0] = 0;
    u[2] = 1;
    u[6] = 1;
    u[8] = 2;

    /*node1,9を除外してGauss EliminationでAu=fを解く*/
    double tmpA[5][5];
    double tmpf[5];

    int list[5] = {1,3,4,5,7};
    for (int i=0; i<5; i++){
        for (int j=0; j<5; j++){
            tmpA[i][j] = largeA[list[i]][list[j]];
        }
        tmpf[i] = largef[list[i]];
    }

    for (int i=0; i<4; i++){
        for (int j=i+1; j<5; j++){
            double aa = tmpA[j][i]/tmpA[i][i];
            tmpf[j] -= aa*tmpf[i];
            for (int k=i+1; k<5; k++) tmpA[j][k] -= aa*tmpA[i][k];
        }
    }
    tmpf[4] /= tmpA[4][4];
    for (int i=3; i>=0; i--){
        for (int j=i+1; j<5; j++) tmpf[i] -= tmpA[i][j]*tmpf[j];
        tmpf[i] /= tmpA[i][i];
    }

    for (int i=0; i<5; i++) u[list[i]] = tmpf[i];

    printf("\n\nAPPROXIMATE VALUES at EACH NODES\n");
    for (int i=0; i<9; i++) printf("u%d: %lf\n", i+1, u[i]);
    return 0;
}