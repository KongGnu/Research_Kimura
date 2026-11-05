#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

//int C,N,O,S,H,X;
typedef struct {
    int C, N, O, S, H, X;
} AtomCount;

//Initial material （他の分子を使いたい場合は以下のコードを書き換えてください）
int Gly = 10;
int Ala = 0;
int Ser = 0;
int Cys = 0;
int X2 = 1;


AtomCount COUNT(int a ,int b, int c, int d, int e){ //炭素原子の行（列）＝0行目〜C-1行目、窒素分子の行（列）＝C行目〜N-1行目、酸素原子の行（列）＝N行目〜O-1行目、硫黄原子の行（列）＝O行目〜S-1行目、水素原子の行（列）＝S行目=H行目、X原子の行（列）＝H+1行目
    //C=0,O=0,H=0,N=0,S=0;
    
    AtomCount counts;

    //Carbon
    counts.C = 2*a + 3*b + 3*c + 3*d;

    //Nitrogen
    counts.N = counts.C + a + b + c + d;

    //Oxygen
    counts.O = counts.N + 2*a + 2*b + 3*c + 2*d;

    //Sulfur
    counts.S = counts.O + d;

    //Hydrogen
    counts.H = counts.S ;

    counts.X = counts.H + 1;

    return counts;
}

void INIT(int a[80][80], AtomCount counts){ //初期分子の行列設定

    //COUNT(Gly,Ala,Ser,Cys,X2);

    //S = O ;
    //O = N ;
    //N = C ;
    //C = 0 ;

    
    int C = 0;
    int N = counts.C;
    int O = counts.N;
    int S = counts.O;
    int H = counts.S;
    int X = counts.X;


    for(int i=0 ; i<=X ; i++){
        for(int j=i ; j<=X ; j++){
            a[i][j]=0;
        }
    }



    //H2NCH2COOH Gly
    for (int i=0 ; i<Gly ; i++ ){
        a[N][H]=2;
        a[C][N]=1;
        a[C][H]=2;
        a[C][C+1]=1;
        a[C+1][O]=2;
        a[C+1][O+1]=1;        
        a[O+1][H]=1;
        C=C+2;
        N=N+1;
        O=O+2;
    }

    //H2NCHCH3COOH Ala
    for (int i=0 ; i<Ala ; i++ ){
        a[N][H]=2;
        a[C][N]=1;
        a[C][H]=1;
        a[C][C+1]=1;
        a[C+1][H]=3;
        a[C][C+2]=1;
        a[C+2][O]=2;
        a[C+2][O+1]=1;
        a[O+1][H]=1;
        C=C+3;
        N=N+1;
        O=O+2;
    }

    //H2NCHCH2OHCOOH Ser
    for (int i=0 ; i<Ser ; i++ ){
        a[N][H]=2;
        a[C][N]=1;
        a[C][H]=1;
        a[C][C+1]=1;
        a[C+1][H]=2;
        a[C+1][O]=1;
        a[O][H]=1;
        a[C][C+2]=1;
        a[C+2][O+1]=2;
        a[C+2][O+2]=1;
        a[O+2][H]=1;      
        C=C+3;
        N=N+1;
        O=O+3;
    }

    //H2NCHCH2SHCOOH Cys
    for (int i=0 ; i<Cys ; i++ ){
        a[N][H]=2;
        a[C][N]=1;
        a[C][H]=1;
        a[C][C+1]=1;
        a[C+1][H]=2;
        a[C+1][S]=1;
        a[S][H]=1;
        a[C][C+2]=1;
        a[C+2][O]=2;
        a[C+2][O+1]=1;
        a[O+1][H]=1;      
        C=C+3;
        N=N+1;
        O=O+2;
        S=S+1;
    }

    //X2
    a[X][X]=X2;

    //H2
    a[H][H]=0;

}

void ENERGY(double b[50][80][80]){ //結合エネルギー

    AtomCount counts = COUNT(Gly,Ala,Ser,Cys,X2);
    int C = counts.C;
    int N = counts.N;
    int O = counts.O;
    int S = counts.S;
    int H = counts.H;
    int X = counts.X;

    for (int i=0 ; i<=X ; i++){
        for (int j=i ;j<=X ; j++){
            b[0][i][j]=0.;
        }
    }

    //C-C結合
    for (int i=0 ; i<C ; i++){
        for (int j=i+1 ; j<C ; j++){

            b[1][i][j]=347.;
            b[2][i][j]=611.;
            b[3][i][j]=837.;
            b[4][i][j]=-1000.;

        }
    }

    //C-N結合
    for (int i=0 ; i<C ; i++){
        for (int j=C ; j<N ; j++){

            b[1][i][j]=305.;
            b[2][i][j]=615.;
            b[3][i][j]=891.;

        }
    }

    //C-O結合
    for (int i=0 ; i<C ;i++){
        for (int j=N ; j<O ; j++){

            b[1][i][j]=358.;
            b[2][i][j]=745.;
            b[3][i][j]=1072.;

        }
    }

    //C-S結合
    for (int i=0 ; i<C ;i++){
        for (int j=O ; j<S ; j++){

            b[1][i][j]=302.;
            b[2][i][j]=482.;
           
        }
    }
    
    //C-H結合
    for (int i=0 ; i<C ; i++){
        b[1][i][H]=414.;
        b[2][i][H]=414.*2.;
        b[3][i][H]=414.*3.;
        b[4][i][H]=414.*4.;
        b[5][i][H]=414.*5.;
    }

    //N-N結合
    for (int i=C ; i<N ; i++){
        for (int j=i+1 ; j<N ; j++){

            b[1][i][j]=161.;
            b[2][i][j]=456.;
            b[3][i][j]=946.;

        }
    }

    //N-O結合
    for (int i=C ; i<N ; i++ ){
        for (int j=N ; j<O ; j++){

            b[1][i][j]=230.;
            b[2][i][j]=598.;

        }
    }

    //N-S結合
    for (int i=C ; i<N ; i++){
        for (int j=O ; j<S ; j++){

            b[1][i][j]=223.;
            b[2][i][j]=444.;
           
        }
    }
    
    //N-H結合
    for (int i=C ; i<N ; i++){
        
        b[1][i][H]=389.;
        b[2][i][H]=389.*2.;
        b[3][i][H]=389.*3.;
        b[4][i][H]=389.*4.;

    }

    //O-O結合
    for (int i=N ; i<O ; i++ ){
        for (int j=i+1 ; j<O ; j++){

            b[1][i][j]=146.;
            b[2][i][j]=498.;

        }
    }

    //O-S結合
    for (int i=N ; i<O ; i++){
        for (int j=O ; j<S ; j++){

            b[1][i][j]=290.;
            b[2][i][j]=594.;
            
        }
    }
    
    //O-H結合
    for (int i=N ; i<O ; i++){
        b[1][i][H]=464.;
        b[2][i][H]=464.*2.;
        b[3][i][H]=464.*3.;

    }

    //S-S結合
    for (int i=O ; i<S ; i++ ){
        for (int j=i+1 ; j<S ; j++){

            b[1][i][j]=230.;
            b[2][i][j]=429.;

        }
    }

    //S-H結合
    for (int i=O ; i<S ; i++){
        b[1][i][H]=374.;
        b[2][i][H]=374.*2.;
        b[3][i][H]=374.*3.;

    }
    
    //H-H結合
    for (int i=1 ; i<50 ; i++){
        b[i][H][H]=436.*(double)i;
    }

    //R-X結合
    for(int i=0 ; i<H ; i++){
        b[1][i][X]=0.;
    }

    //H-X結合
    //X-X結合
    for(int i=1 ; i<50 ; i++){
        b[i][H][X]=0.;
        b[i][X][X]=-1000000.*(double)i;
    }
}

void LENGTH(double b[4][80][80]){ //結合長 [pm]

    AtomCount counts = COUNT(Gly,Ala,Ser,Cys,X2);
    int C = counts.C;
    int N = counts.N;
    int O = counts.O;
    int S = counts.S;
    int H = counts.H;
    int X = counts.X;

    for (int i=0 ; i<=X ; i++){
        for (int j=i ;j<=X ; j++){
            b[0][i][j]=0.;
        }
    }

    //C-C結合
    for (int i=0 ; i<C ; i++){
        for (int j=i+1 ; j<C ; j++){

            b[1][i][j]=75.+75.;
            b[2][i][j]=67.+67.;
            b[3][i][j]=60.+60.;

        }
    }

    //C-N結合
    for (int i=0 ; i<C ; i++){
        for (int j=C ; j<N ; j++){

            b[1][i][j]=75.+71.;
            b[2][i][j]=67.+60.;
            b[3][i][j]=60.+54.;

        }
    }

    //C-O結合
    for (int i=0 ; i<C ;i++){
        for (int j=N ; j<O ; j++){

            b[1][i][j]=75.+63.;
            b[2][i][j]=67.+57.;
            b[3][i][j]=60.+53.;

        }
    }

    //C-S結合
    for (int i=0 ; i<C ;i++){
        for (int j=O ; j<S ; j++){

            b[1][i][j]=75.+103.;
            b[2][i][j]=67.+94.;
           
        }
    }
    
    //C-H結合
    for (int i=0 ; i<C ; i++){
        b[1][i][H]=75.+32.;
    }

    //C-X結合
    for (int i=0 ; i<C ; i++){
        b[1][i][X]=75.;
    }

    //N-N結合
    for (int i=C ; i<N ; i++){
        for (int j=i+1 ; j<N ; j++){

            b[1][i][j]=71.+71.;
            b[2][i][j]=60.+60.;
            b[3][i][j]=54.+54.;

        }
    }

    //N-O結合
    for (int i=C ; i<N ; i++ ){
        for (int j=N ; j<O ; j++){

            b[1][i][j]=71.+63.;
            b[2][i][j]=60.+57.;

        }
    }

    //N-S結合
    for (int i=C ; i<N ; i++){
        for (int j=O ; j<S ; j++){

            b[1][i][j]=71.+103.;
            b[2][i][j]=60.+94.;
           
        }
    }
    
    //N-H結合
    for (int i=C ; i<N ; i++){
        
        b[1][i][H]=71.+32.;

    }

    //N-X結合
    for (int i=C ; i<N ; i++){
        b[1][i][X]=71.;
    }

    //O-O結合
    for (int i=N ; i<O ; i++ ){
        for (int j=i+1 ; j<O ; j++){

            b[1][i][j]=63.+63.;
            b[2][i][j]=57.+57.;

        }
    }

    //O-S結合
    for (int i=N ; i<O ; i++){
        for (int j=O ; j<S ; j++){

            b[1][i][j]=63.+103.;
            b[2][i][j]=57.+94.;
            
        }
    }
    
    //O-H結合
    for (int i=N ; i<O ; i++){
        b[1][i][H]=63.+32.;

    }

    //O-X結合
    for (int i=N ; i<O ; i++){
        b[1][i][X]=63.;
    }

    //S-S結合
    for (int i=O ; i<S ; i++ ){
        for (int j=i+1 ; j<S ; j++){

            b[1][i][j]=103.+103.;
            b[2][i][j]=94.+94.;

        }
    }

    //S-H結合
    for (int i=O ; i<S ; i++){
        b[1][i][H]=103.+32.;

    }

    //S-X結合
    for (int i=O ; i<S ; i++){
        b[1][i][X]=103.;
    }
    
    //H-H結合
    b[1][H][H]=32.+32.;

    //H-X結合
    b[1][H][X]=32.;

    //X-X結合は0

}

struct nonzero{ 
    int row1;
    int col1;
    int row2;
    int col2;
    int row3;
    int col3;
    int row4;
    int col4;
};

// Residue structure for peptide detection
struct Residue{
    int carboxyC; // carbon index of carboxyl C (= C=O)
    int alphaC;   // alpha carbon index
    int aminoN;   // nitrogen index of amino group
};

typedef struct{
    int disk[80];
    int low[80];
    int parent[80];
    int visited[80]; //DFS用の訪問済みフラグ

    int timer;

    int bridge[80][80]; //DFSで見つかった橋の情報を格納する配列
}SYSTEM;

typedef struct{
    double mass;
    int molecule;
}ATOM;

typedef struct{
    int id;
    int nAtom;
    int nBond;
    double mw;
    double theta;  //振動特性温度　
}MOLECULE;

//ATOM atom[80];
//MOLECULE molecule[80];

void MASS(ATOM atom[80]){ //原子の質量と分子の分子量を設定

    AtomCount counts = COUNT(Gly,Ala,Ser,Cys,X2);
    int C = counts.C;
    int N = counts.N;
    int O = counts.O;
    int S = counts.S;
    int H = counts.H;
    int X = counts.X;

    for(int i=0 ; i<=X ; i++){
        atom[i].molecule=0;
        atom[i].mass=0.;
    }

    for(int i=0 ; i<C ; i++){
        atom[i].mass=12.011; //C
    }
    for(int i=C ; i<N ; i++){
        atom[i].mass=14.007; //N
    }
    for(int i=N ; i<O ; i++){
        atom[i].mass=15.999; //O
    }
    for(int i=O ; i<S ; i++){
        atom[i].mass=32.06;  //S
    }
    for(int i=H ; i<H+1 ; i++){
        atom[i].mass=1.008;  //H
    }

    //Xの質量は0とする

}

void DFS(SYSTEM *sys, ATOM atom[80], MOLECULE molecule[80], int now, int mol, int b[80][80], int H){
    int next;
    sys->visited[now]=1;

    atom[now].molecule=mol;
    molecule[mol].mw=molecule[mol].mw+atom[now].mass;
    molecule[mol].nAtom++;

    sys->disk[now]=sys->low[now]=++sys->timer;

    for(next=0; next<=H; next++){
        if(now<next && b[now][next]==0) continue; //隣接していない場合はスキップ
        else if(next<now && b[next][now]==0) continue; //隣接していない場合はスキップ

        if(sys->parent[now]==next) continue; //親への辺は無視

        if(sys->visited[next]==0){
            sys->parent[next]=now;
            DFS(sys, atom, molecule, next, mol, b, H);

            if(sys->low[next]<sys->low[now]){
                sys->low[now]=sys->low[next];
            }

            //橋の判定
            if(now<next && b[now][next]==1){
                if(sys->low[next]>sys->disk[now]){
                    sys->bridge[now][next]=1;
                }
            }
            else if(next<now && b[next][now]==1){
                if(sys->low[next]>sys->disk[now]){
                    sys->bridge[next][now]=1;
                }
            }

        }else{
            if(sys->disk[next]<sys->low[now]){
                sys->low[now]=sys->disk[next];
            }
        }
    }
}

void BOND(ATOM atom[80], MOLECULE molecule[80], int mol, int b[80][80], int C, int N, int O, int S, int H){
    for(int i=0; i<=H; i++){
        for(int j=i+1; j<=H; j++){
            if((atom[i].molecule == mol) && (atom[j].molecule == mol)){
                if(b[i][j]>0){

                    if(i<C){
                        if(j<C){
                            if(b[i][j]==1){ //C-C　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1403.;
                                continue;
                            }else if(b[i][j]==2){ //C=C　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 2374.;
                                continue;
                            }else if(b[i][j]==3){ //C三C　三重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 3137.;
                                continue;
                            }
                        }else if(j<N){
                            if(b[i][j]==1){ //C-N　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1669.;
                                continue;
                            }else if(b[i][j]==2){ //C=N　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 2319.;
                                continue;
                            }else if(b[i][j]==3){ //C三N　三重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 3230.;
                                continue;
                            }
                        }else if(j<O){
                            if(b[i][j]==1){ //C-O　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1760.;
                                continue;
                            }else if(b[i][j]==2){ //C=O　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 2460.;
                                continue;
                            }else if(b[i][j]==3){ //C三O　三重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 3122.;
                                continue;
                            }
                        }else if(j<S){
                            if(b[i][j]==1){ //C-S　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 960.;
                                continue;
                            }else if(b[i][j]==2){ //C=S　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1592.;
                                continue;
                            }
                        }else{  //C-H
                            molecule[mol].nBond = molecule[mol].nBond + b[i][H];
                            molecule[mol].theta = molecule[mol].theta + b[i][H]*4180.;
                            continue;
                        }

                    }else if(i<N){
                        if(j<N){
                            if(b[i][j]==1){ //N-N　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1511.;
                                continue;
                            }else if(b[i][j]==2){ //N=N　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 2306.;
                                continue;
                            }else if(b[i][j]==3){ //N三N　三重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 3393.;
                                continue;
                            }
                        }else if(j<O){
                            if(b[i][j]==1){ //N-O　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1378.;
                                continue;
                            }else if(b[i][j]==2){ //N=O　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 2269.;
                                continue;
                            }
                        }else if(j<S){
                            if(b[i][j]==1){ //N-S　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1432.;
                                continue;
                            }else if(b[i][j]==2){ //N=S　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1732.;
                                continue;
                            }
                        }else{  //N-H
                            molecule[mol].nBond = molecule[mol].nBond + b[i][H];
                            molecule[mol].theta = molecule[mol].theta + b[i][H]*4820.;
                            continue;
                        }

                    }else if(i<O){
                        if(j<O){
                            if(b[i][j]==1){ //O-O　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1262.;
                                continue;
                            }else if(b[i][j]==2){ //O=O　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 2274.;
                                continue;
                            }
                        }else if(j<S){
                            if(b[i][j]==1){ //O-S　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 935.;
                                continue;
                            }else if(b[i][j]==2){ //O=S　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1619.;
                                continue;
                            }
                        }else{  //O-H
                            molecule[mol].nBond = molecule[mol].nBond + b[i][H];
                            molecule[mol].theta = molecule[mol].theta + b[i][H]*5208.;
                            continue;
                        }

                    }else if(i<S){
                        if(j<S){
                            if(b[i][j]==1){ //S-S　単結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 764.;
                                continue;
                            }else if(b[i][j]==2){ //S=S　二重結合
                                molecule[mol].nBond++;
                                molecule[mol].theta = molecule[mol].theta + 1007.;
                                continue;
                            }
                        }else{  //S-H
                            molecule[mol].nBond = molecule[mol].nBond + b[i][H];
                            molecule[mol].theta = molecule[mol].theta + b[i][H]*3957.;
                            continue;
                        }

                    }else{  //H-H
                        molecule[mol].nBond = molecule[mol].nBond + b[H][H];
                        molecule[mol].theta = molecule[mol].theta + b[H][H]*6331.;
                    }
                }
            }
        }
    }
    molecule[mol].theta = molecule[mol].theta/molecule[mol].nBond;
}

int max(int a, int b){
    if(a<b){
        return a;
    }else{
        return b;
    }
}


int main() {
    FILE* f = fopen("file_abundance_10Gly_test_calculation_solid_entropy.txt", "w");//アバンダンス出力ファイル
    FILE* s = fopen("file_marvinview_10Gly_test_calculation_solid_entropy.sdf", "w");//Marvin View用のSDFファイル

    //Parameters
    int TRIAL=100000;    //試行回数
    double T=300.;      //温度 K
    int STEP=240000;      //全ステップ数
    int PostUV=120000;    //PostUVphase開始ステップ
    double alpha=0.5;
    double beta=300.; //[kJ/mol]
    //double timescale=0.000019*3e-2/600; //yr
    //double Ea=40*(T/100)*(1.1+0.048*log10((timescale/1000.)*(T/100)));//[kJ/mol]：ある温度Tにおいて、timescaleの期間内に起きる反応の活性化エネルギーの最大値
    //if(Ea<0) Ea=0.;//[kJ/mol]：無反応の活性化エネルギーは0とする

    //分子種判別　カウント間隔
    int interval=100;

    //衝突判定パラメーター
    double dt=3.0e-2;          //衝突確率を求めるための時間刻み [s]
    double dt_post=3.0e-2;    //PostUVphaseの時間刻み [s]

    double f_ch=3.28; //衝突頻度 [1/s]
    double prob_ch=1.0-exp(-f_ch*dt); //荷電粒子の衝突確率
    //prob_ch=0; //荷電粒子の衝突確率を固定

    double f_th=3.37; //衝突頻度 [1/s]
    double prob_th=1.0-exp(-f_th*dt); //衝突確率
    double prob_th_post=1.0-exp(-f_th*dt_post); //PostUVphaseの衝突確率
    //double prob_th=1.; //衝突確率を固定　毎ステップ何かしらの反応が起こるとする
    //double prob_th_post=1.; //PostUVphaseの衝突確率を固定

    static double B_energy[50][80][80]; //結合エネルギー
    static double B_length[4][80][80]; //結合長
    ENERGY(B_energy);
    LENGTH(B_length);

    /*int B_matrix_init[80][80]; //初期結合行列
    INIT(B_matrix_init);

    //Ramdom number
    double r0=7788776655. ;
    double ra=48271. ; 
    double rm=2147483647. ;*/

    //気体定数
    double R=0.00831; 

    //Abundance
    static int aacid_abundance[30000];
    //static double sugar_abundance[18000];
    //static double CO_abundance[18000];
    //static double CO2_abundance[18000];
    //static double H2_abundance[18000];
    //static double H2O_abundance[18000];
    //static double peptidelike_abundance[18000];
    static int peptide_dimer_abundance[30000];
    static int peptide_trimer_abundance[30000];
    static int peptide_oligomer_abundance[30000];

    static int collision_ch[30000]; //荷電粒子の衝突の回数
    static int collision_th[30000]; //衝突の回数

    int aacid_abundance_UV=0;
    int peptide_dimer_abundance_UV=0;
    int peptide_trimer_abundance_UV=0;
    int peptide_oligomer_abundance_UV=0;
    int collision_ch_UV=0;
    int collision_th_UV=0;
    
    int aacid_abundance_PostUV=0;
    int peptide_dimer_abundance_PostUV=0;
    int peptide_trimer_abundance_PostUV=0;
    int peptide_oligomer_abundance_PostUV=0;
    int collision_ch_PostUV=0;
    int collision_th_PostUV=0;

    for (int i=0 ; i<30000 ; i++){
        aacid_abundance[i]=0;
        //sugar_abundance[i]=0;
        //H2_abundance[i]=0;
        //H2O_abundance[i]=0;
        //CO_abundance[i]=0;
        //CO2_abundance[i]=0;
        //peptidelike_abundance[i]=0;
        peptide_dimer_abundance[i]=0;
        peptide_trimer_abundance[i]=0;
        peptide_oligomer_abundance[i]=0;
        collision_ch[i]=0;
        collision_th[i]=0;
    }

    //MarvinViewファイル用の変数
    int abc=0;

    // ★ OpenMP並列化開始
    #pragma omp parallel reduction(+:aacid_abundance[:30000], peptide_dimer_abundance[:30000], peptide_trimer_abundance[:30000], peptide_oligomer_abundance[:30000], collision_ch[:30000], collision_th[:30000], aacid_abundance_UV, peptide_dimer_abundance_UV, peptide_trimer_abundance_UV, peptide_oligomer_abundance_UV, collision_ch_UV, collision_th_UV, aacid_abundance_PostUV, peptide_dimer_abundance_PostUV, peptide_trimer_abundance_PostUV, peptide_oligomer_abundance_PostUV, collision_ch_PostUV, collision_th_PostUV) private(abc)
    {
        AtomCount counts = COUNT(Gly,Ala,Ser,Cys,X2);

        int C = counts.C;
        int N = counts.N;
        int O = counts.O;
        int S = counts.S;
        int H = counts.H;
        int X = counts.X;

        #pragma omp for
        for (int I=0 ; I <TRIAL; I++){//試行回数

            // ★ 大きな配列を動的割り当て（ヒープに確保）
            enum { B_NZ_CAP = 10000, R_NZ_CAP = 100000 };
#define CHECK_CAPACITY(index, cap, name) do { \
                if ((index) >= (cap)) { \
                    fprintf(stderr, "%s overflow in trial %d\n", (name), I); \
                    free(B_matrix_local); \
                    free(B_nz_local); \
                    free(R_nz_local); \
                    free(D_energy_local); \
                    free(dH_local); \
                    free(dS_local); \
                    free(P_local); \
                    exit(EXIT_FAILURE); \
                } \
            } while (0)
            int (*B_matrix_local)[80] = malloc(80 * 80 * sizeof(*B_matrix_local));
            struct nonzero *B_nz_local = malloc(B_NZ_CAP * sizeof(*B_nz_local));
            struct nonzero *R_nz_local = malloc(R_NZ_CAP * sizeof(*R_nz_local));
            double *D_energy_local = malloc(R_NZ_CAP * sizeof(*D_energy_local));
            double *dH_local = malloc(R_NZ_CAP * sizeof(*dH_local));
            double *dS_local = malloc(R_NZ_CAP * sizeof(*dS_local));
            double *P_local = malloc(R_NZ_CAP * sizeof(*P_local));

            if (B_matrix_local == NULL || B_nz_local == NULL || R_nz_local == NULL || D_energy_local == NULL || dH_local == NULL || dS_local == NULL || P_local == NULL) {
                fprintf(stderr, "malloc failed in trial %d\n", I);
                free(B_matrix_local);
                free(B_nz_local);
                free(R_nz_local);
                free(D_energy_local);
                free(dH_local);
                free(dS_local);
                free(P_local);
                exit(EXIT_FAILURE);
            }

            for(int i=0 ; i<R_NZ_CAP ; i++){
                D_energy_local[i]=0.;
                dH_local[i]=0.;
                dS_local[i]=0.;
                P_local[i]=0.;
            }
            
            double ra=48271.;
            double rm=2147483647.;
            double r0_local = 7788776655.0 + I * 12345.0;  // スレッド個別初期値

            // ★ INIT呼び出しを #pragma omp critical で保護（グローバル変数C,N,O,S,H,Xへのアクセスを同期化）
            //#pragma omp critical
            //{
                INIT(B_matrix_local, counts); //結合行列を初期物質にリセット
            //}

            for (int J=0 ; J<STEP; J++){ //ステップ回数

                if(J<PostUV){ 
                    B_matrix_local[X][X]=X2; //UVphaseの間はX2を含める
                }else{
                    B_matrix_local[X][X]=0; //PostUVphaseはX2なし
                }

                ATOM atom[80];
                ATOM test_atom[80];
                MOLECULE molecule[80];
                MOLECULE test_mol[80];

                MASS(atom);
                MASS(test_atom);

                for(int i=0 ; i<H ; i++){
                    molecule[i].id=i;
                    molecule[i].nAtom=0;
                    molecule[i].mw=0.;
                    test_mol[i].id=i;
                    test_mol[i].nAtom=0;
                    test_mol[i].mw=0.;
                }

                //static struct nonzero B_nz[10000]; //ノンゼロ(行列成分のうち、結合を表す成分)リスト

                struct nonzero CO_dou[100]; //CO二重結合のリスト
                struct nonzero CO_tri[100]; //CO三重結合のリスト
                struct nonzero SO2_mol[100]; //SO二重結合のリスト
                struct nonzero ROR_mol[100]; //R-O-R'構造のリスト

                int Total_nz=0; //ノンゼロ成分の個数
                int CO2=0, CO3=0; //CO二重結合, CO三重結合の個数
                int SO2=0; //SO二重結合の個数
                int ROR=0; //R-O-R'構造の個数

                int bb=0; //R行列の数＝反応候補の数

                int S_bond[80]; //硫黄原子の結合手の数 []内は行列サイズ

                for(int i=0; i<=X; i++){
                    S_bond[i]=0;    //硫黄原子以外については結合数0が記録される
                }

                for(int i=0; i<=X; i++){
                    for(int j=O; j<S; j++){
                        if(i<j){
                            S_bond[j]=S_bond[j]+B_matrix_local[i][j];
                        }
                        if(j<i){
                            S_bond[j]=S_bond[j]+B_matrix_local[j][i];
                        }
                    }
                }

                //--------------------------------荷電粒子による反応候補の選択--------------------------------
                r0_local = fmod(ra*r0_local,rm); //乱数

                if(B_matrix_local[X][X]>=1 && r0_local/rm<=prob_ch){

                    //ノンゼロ成分の行(row1)、列(col1)を記録
                    for(int i=0;i<=H;i=i+1){ //C~H列の探索(X2はリストに含めない)
                        for(int j=i; j<=X; j=j+1){ //C~X行の探索
                            if(S_bond[i]>2 || S_bond[j]>2){ //S原子の結合手が3本以上のときは、反応候補から除外する
                                continue;
                            }
                            if(j==H || j==X){ 
                                for (int k=0 ; k<B_matrix_local[i][j] ; k++){ 
                                    CHECK_CAPACITY(Total_nz, B_NZ_CAP, "B_nz_local");
                                    B_nz_local[Total_nz].row1=i;
                                    B_nz_local[Total_nz].col1=j;
                                    Total_nz=Total_nz+1;
                                }
                            }else if (i<C && j>=N){ 
                                if(B_matrix_local[i][j]==1 || B_matrix_local[i][j]==2){//CO分子はカウントしないように
                                    CHECK_CAPACITY(Total_nz, B_NZ_CAP, "B_nz_local");
                                    B_nz_local[Total_nz].row1=i;
                                    B_nz_local[Total_nz].col1=j;
                                    Total_nz=Total_nz+1;
                                }
                            }else if(B_matrix_local[i][j]!=0){
                                CHECK_CAPACITY(Total_nz, B_NZ_CAP, "B_nz_local");
                                B_nz_local[Total_nz].row1=i;
                                B_nz_local[Total_nz].col1=j;
                                Total_nz=Total_nz+1;
                            }
                        }
                    }

                    //CO二重結合、三重結合成分の行(row1)、列(col1)を記録
                    for(int i=0;i<C;i=i+1){
                        for(int j=N; j<O; j=j+1){
                            if (B_matrix_local[i][j]==2){
                                CO_dou[CO2].row1=i;
                                CO_dou[CO2].col1=j;
                                CO2=CO2+1;
                            }else if (B_matrix_local[i][j]==3){
                                CO_tri[CO3].row1=i;
                                CO_tri[CO3].col1=j;
                                CO3=CO3+1;
                            }
                        }
                    }

                    //SO二重結合の行(row1)、列(col1)を記録
                    for(int i=N;i<O;i=i+1){
                        for(int k=O; k<S; k=k+1){
                            if (B_matrix_local[i][k]==2){
                                SO2_mol[SO2].row1=i;    //O原子の番号
                                SO2_mol[SO2].col1=k;    //S原子の番号
                                SO2=SO2+1;
                            }
                        }
                    }

                    //Sの結合数が4または6で、結合数が変わらない反応の行列成分を記録
                        for(int i=0; i<H; i++){
                            int count=0;
                            if(S_bond[i]==4){
                                for(int l=0; l<SO2; l++){
                                    if (SO2_mol[l].col1==i){   //O=S(-R)-R' 
                                        if(SO2_mol[l].col1!=count){//ノンゼロ成分の多重カウントを防止
                                            for(int j=0; j<=X; j++){
                                                if(j==H || j==X){
                                                    for (int k=0 ; k<B_matrix_local[SO2_mol[l].col1][j] ; k++){
                                                        B_nz_local[Total_nz].row1=SO2_mol[l].col1;    //S原子の番号
                                                        B_nz_local[Total_nz].col1=j;                  //HまたはX原子の番号
                                                        Total_nz=Total_nz+1;
                                                        count=i;
                                                    }
                                                }else if(SO2_mol[l].col1<j){
                                                    if(B_matrix_local[SO2_mol[l].col1][j]!=0){
                                                        B_nz_local[Total_nz].row1=SO2_mol[l].col1;    //S原子の番号
                                                        B_nz_local[Total_nz].col1=j;                  //RまたはR'原子の番号
                                                        Total_nz=Total_nz+1;
                                                        count=i;
                                                    }
                                                }else if(j<SO2_mol[l].col1){
                                                    if(B_matrix_local[j][SO2_mol[l].col1]==1){
                                                        B_nz_local[Total_nz].row1=j;                  //RまたはR'原子の番号
                                                        B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                        Total_nz=Total_nz+1;
                                                        count=i;
                                                    }else if (B_matrix_local[j][SO2_mol[l].col1]==2){
                                                        if(j<N && O<=j){   //jはO原子を除く
                                                            B_nz_local[Total_nz].row1=j;                  //RまたはR'原子の番号
                                                            B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                            Total_nz=Total_nz+1;
                                                            count=i;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    for(int m=l+1; m<SO2; m++){
                                        if(SO2_mol[l].col1==i && SO2_mol[m].col1==i){   //O=S=O'

                                            B_nz_local[Total_nz].row1=SO2_mol[l].row1;    //O原子の番号
                                            B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                            Total_nz=Total_nz+1;

                                            B_nz_local[Total_nz].row1=SO2_mol[m].row1;    //O'原子の番号
                                            B_nz_local[Total_nz].col1=SO2_mol[m].col1;    //S原子の番号
                                            Total_nz=Total_nz+1;
                                            //break;
                                        }
                                    }
                                }
                            }else if (S_bond[i]==6){
                                for(int l=0; l<SO2; l++){
                                    for(int m=l+1; m<SO2; m++){
                                        if (SO2_mol[l].col1==i && SO2_mol[m].col1==i){   //O=S(-R)(=O')-R'
                                            if(SO2_mol[l].col1!=count){//ノンゼロ成分の多重カウントを防止
                                                for(int j=0; j<=X; j++){
                                                    if(j==H || j==X){
                                                        for (int k=0 ; k<B_matrix_local[SO2_mol[l].col1][j] ; k++){
                                                            B_nz_local[Total_nz].row1=SO2_mol[l].col1;    //S原子の番号
                                                            B_nz_local[Total_nz].col1=j;                  //HまたはX原子の番号
                                                            Total_nz=Total_nz+1;
                                                            count=i;
                                                        }
                                                    }else if(SO2_mol[l].col1<j){
                                                        if(B_matrix_local[SO2_mol[l].col1][j]!=0){
                                                            if(S_bond[j]==6){   //Sの結合数が6のときのみ反応候補に加える
                                                                B_nz_local[Total_nz].row1=SO2_mol[l].col1;    //S原子の番号
                                                                B_nz_local[Total_nz].col1=j;                  //RまたはR'原子の番号
                                                                Total_nz=Total_nz+1;
                                                            }
                                                            count=i;
                                                        }
                                                    }else if(j<SO2_mol[l].col1 && j<O){
                                                        if(B_matrix_local[j][SO2_mol[l].col1]==1){
                                                            B_nz_local[Total_nz].row1=j;                  //RまたはR'原子の番号
                                                            B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                            Total_nz=Total_nz+1;
                                                            count=i;
                                                        }else if (B_matrix_local[j][SO2_mol[l].col1]==2){
                                                            if(j<N && O<=j){   //jはO原子を除く
                                                                B_nz_local[Total_nz].row1=j;                  //RまたはR'原子の番号
                                                                B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                                Total_nz=Total_nz+1;
                                                                count=i;
                                                            }
                                                        }
                                                    }else{
                                                        count=i;
                                                    }
                                                }
                                            }
                                        }
                                        for(int n=m+1; n<SO2; n++){
                                            if(SO2_mol[l].col1==i && SO2_mol[m].col1==i && SO2_mol[n].col1==i){   //O=S(=O')=O"

                                                B_nz_local[Total_nz].row1=SO2_mol[l].row1;    //O原子の番号
                                                B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                Total_nz=Total_nz+1;

                                                B_nz_local[Total_nz].row1=SO2_mol[m].row1;    //O'原子の番号
                                                B_nz_local[Total_nz].col1=SO2_mol[m].col1;    //S原子の番号
                                                Total_nz=Total_nz+1;

                                                B_nz_local[Total_nz].row1=SO2_mol[n].row1;    //O"原子の番号
                                                B_nz_local[Total_nz].col1=SO2_mol[n].col1;    //S原子の番号
                                                Total_nz=Total_nz+1;
                                                //break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    //*/

                    /*/R-O-R'構造の行、列を記録
                    for(int i=N; i<O; i=i+1){
                        for(int j=0; j<=H; j=j+1){
                            for(int k=j+1; k<=X; k=k+1){
                                if(i<j){
                                    if(B_matrix_local[i][j]==1 && B_matrix_local[i][k]==1){
                                        ROR_mol[ROR].row1=i;    //O原子の番号
                                        ROR_mol[ROR].col1=j;    //Rの番号
                                        ROR_mol[ROR].row2=i;    //O原子の番号
                                        ROR_mol[ROR].col2=k;    //R'の番号
                                        ROR=ROR+1;
                                    }
                                }
                                if(j<i && i<k){
                                    if(B_matrix_local[j][i]==1 && B_matrix_local[i][k]==1){
                                        ROR_mol[ROR].row1=i;    //O原子の番号
                                        ROR_mol[ROR].col1=j;    //Rの番号
                                        ROR_mol[ROR].row2=i;    //O原子の番号
                                        ROR_mol[ROR].col2=k;    //R'の番号
                                        ROR=ROR+1;
                                    }
                                }
                                if(k<i){
                                    if(B_matrix_local[j][i]==1 && B_matrix_local[k][i]==1){
                                        ROR_mol[ROR].row1=i;    //O原子の番号
                                        ROR_mol[ROR].col1=j;    //Rの番号
                                        ROR_mol[ROR].row2=i;    //O原子の番号
                                        ROR_mol[ROR].col2=k;    //R'の番号
                                        ROR=ROR+1;
                                    }
                                }
                            }
                        }

                        if(B_matrix_local[i][H]==2){
                            ROR_mol[ROR].row1=i;    //O原子の番号
                            ROR_mol[ROR].col1=H;    //H原子の番号
                            ROR_mol[ROR].row2=i;    //O原子の番号
                            ROR_mol[ROR].col2=H;    //H原子の番号
                            ROR=ROR+1;
                        }
                        if(B_matrix_local[i][X]==2){
                            ROR_mol[ROR].row1=i;    //O原子の番号
                            ROR_mol[ROR].col1=X;    //X原子の番号
                            ROR_mol[ROR].row2=i;    //O原子の番号
                            ROR_mol[ROR].col2=X;    //X原子の番号
                            ROR=ROR+1;
                        }
                    }//*/


                    //static struct nonzero R_nz[100000]; //R行列のノンゼロリスト（切断する結合＝-1、生成する結合＝+1）
                    //int bb=0; //R行列の数＝反応候補の数
                    //R_nz[bb].row1, R_nz[bb].col1: 切断する結合の行、列番号
                    //R_nz[bb].row2, R_nz[bb].col2: 切断するもう一つの結合の行、列番号
                    //R_nz[bb].row3, R_nz[bb].col3: 生成する結合の行、列番号
                    //R_nz[bb].row4, R_nz[bb].col4: 生成するもう一つの結合の行、列番号

                    //int ch=0; //荷電粒子による反応候補の数

                    //static double D_energy[100000];//各反応候補の活性化エネルギー
                    //static double dH[100000];//各反応候補のエンタルピー変化

                    if(J==PostUV-1){
                        collision_ch_UV++; //衝突の回数をカウント
                    }else if((STEP!=PostUV) && (J==STEP-1)){
                        collision_ch_PostUV++; //衝突の回数をカウント
                    }else if(J%interval==0){
                        collision_ch[J/interval]=collision_ch[J/interval]+1; //衝突の回数をカウント
                    }

                    //光解離反応 R-R' + X2 -> R-X + R'-X
                    for (int l=0 ; l<Total_nz ; l++){
                        for (int j=0 ; j<B_matrix_local[X][X] ; j++){
                            if(B_nz_local[l].col1!=X){
                                CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1;
                                R_nz_local[bb].row2=X,            R_nz_local[bb].col2=X;
                                R_nz_local[bb].row3=B_nz_local[l].row1, R_nz_local[bb].col3=X;
                                R_nz_local[bb].row4=B_nz_local[l].col1, R_nz_local[bb].col4=X;
                                D_energy_local[bb]=0;//活性化エネルギーは0と仮定
                                dH_local[bb]=100.;//100は適当な値です
                                bb=bb+1;
                            }
                        }
                    }

                    //CO分解反応 CO + X-X -> X-CO-X
                    for(int i=0 ; i<CO3 ; i++){
                        for (int k=0 ; k<B_matrix_local[X][X] ; k++){
                            CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                            R_nz_local[bb].row1=X             , R_nz_local[bb].col1=X;
                            R_nz_local[bb].row2=CO_tri[i].row1, R_nz_local[bb].col2=CO_tri[i].col1;
                            R_nz_local[bb].row3=CO_tri[i].row1, R_nz_local[bb].col3=X;
                            R_nz_local[bb].row4=CO_tri[i].row1, R_nz_local[bb].col4=X;
                            D_energy_local[bb]=0;//活性化エネルギーは0と仮定
                            dH_local[bb]=100.;//100は適当な値です
                            bb = bb + 1;
                        }
                    }

                    int ch_so=bb;

                    //SO2分解反応 SO2 + X-X -> X-O-X + SO
                    for(int m=0; m<=X; m++){
                        if(S_bond[m]==4 || S_bond[m]==6){ //Sの結合数が4または6のとき
                            for(int i=0; i<SO2; i++){
                                if(SO2_mol[i].col1==m){
                                    for (int k=0 ; k<B_matrix_local[X][X] ; k++){
                                        CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                        R_nz_local[bb].row1=X              , R_nz_local[bb].col1=X;
                                        R_nz_local[bb].row2=SO2_mol[i].row1, R_nz_local[bb].col2=SO2_mol[i].col1;
                                        R_nz_local[bb].row3=SO2_mol[i].row1, R_nz_local[bb].col3=X;
                                        R_nz_local[bb].row4=SO2_mol[i].row1, R_nz_local[bb].col4=X;
                                        D_energy_local[bb]=0;//活性化エネルギーは0と仮定
                                        dH_local[bb]=100.;//100は適当な値です
                                        bb = bb + 1;
                                    }
                                }
                            }
                        }
                    }

                    //行列の書き換え

                    int ch = bb; //荷電粒子による反応候補の数

                    //反応速度の和
                    double Total_exp_ch = 0;
                    for (int i=0 ; i<ch ; i++){
                        Total_exp_ch = Total_exp_ch + exp(-(D_energy_local[i])/(R*T));
                    }

                    //ランダムに反応選択 → E行列（生成物行列）生成
                    int re_ch=0; //選ばれた反応の番号
                    double ex_ch = 0 ;
                    r0_local = fmod(ra*r0_local,rm); //乱数
                    for (int i=0 ; i<ch ; i++ ){
                        if (ex_ch <= r0_local/rm && r0_local/rm < ex_ch + (exp(-(D_energy_local[i])/(R*T))/Total_exp_ch)){
                            re_ch=i;//選ばれた反応の番号
                            B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] + 1;//made
                            B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] + 1;//made
                            B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] = B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] - 1;//broken
                            B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] = B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] - 1;//broken
                            break;
                        }
                        ex_ch = ex_ch + exp(-(D_energy_local[i])/(R*T))/Total_exp_ch;
                    }
                }

                //--------------------------------荷電粒子による反応候補の選択終了--------------------------------*/

                for(int i=0 ; i<R_NZ_CAP ; i++){
                    D_energy_local[i]=0.;
                    dH_local[i]=0.;
                    dS_local[i]=0.;
                    P_local[i]=0.;
                }

                SYSTEM sys; //分子の情報を更新
                SYSTEM test_sys; //分子の情報を更新

                sys.timer=0;
                test_sys.timer=0;

                for(int i=0; i<H; i++){
                    sys.disk[i]=0;
                    sys.low[i]=0;
                    sys.parent[i]=0;
                    sys.visited[i]=0; //分子の探索済みフラグの初期化

                    test_sys.disk[i]=0;
                    test_sys.low[i]=0;
                    test_sys.parent[i]=0;
                    test_sys.visited[i]=0; //分子の探索済みフラグの初期化

                    for(int j=0; j<H; j++){ //橋の有無の初期化
                        sys.bridge[i][j]=0;
                        test_sys.bridge[i][j]=0;
                    }
                }

                int current_molecule=0; //現在の分子番号

                for(int i=0; i<H; i++){
                    if(sys.visited[i]==1) continue;

                    molecule[current_molecule].id=current_molecule; //分子番号の記録
                    molecule[current_molecule].mw=0.0; //分子量の初期化
                    molecule[current_molecule].nAtom=0; //分子内の原子数の初期化
                    molecule[current_molecule].nBond=0; //分子内の結合数の初期化
                    molecule[current_molecule].theta=0; //分子内の結合の振動特性温度の初期化

                    DFS(&sys, atom, molecule, i, current_molecule, B_matrix_local, H); //深さ優先探索で分子の情報を更新

                    current_molecule++;
                }

                for(int i=0; i<current_molecule; i++){
                    BOND(atom, molecule, i, B_matrix_local, C, N, O, S, H);
                }

                //--------------------------------衝突による反応候補の選択--------------------------------
                r0_local = fmod(ra*r0_local,rm); //乱数

                if((J<PostUV && r0_local/rm<=prob_th) || (J>=PostUV && r0_local/rm<=prob_th_post)){   //衝突が起きた場合にのみ、反応候補リストを作成する

                    /*if(J%interval==0){
                        collision_th[J/interval]=collision_th[J/interval]+1; //衝突の回数をカウント
                    }*/


                    //結合行列のノンゼロ成分を再度カウント

                    Total_nz=0; //ノンゼロ成分の個数
                    CO2=0, CO3=0; //CO二重結合, CO三重結合の個数
                    SO2=0; //SO二重結合の個数
                    ROR=0; //R-O-R'構造の個数

                    bb=0; //R行列の数＝反応候補の数

                    for(int i=0; i<=X; i++){
                        S_bond[i]=0;    //硫黄原子以外については結合数0が記録される
                    }

                    for(int i=0; i<=X; i++){
                        for(int j=O; j<S; j++){
                            if(i<j){
                                S_bond[j]=S_bond[j]+B_matrix_local[i][j];
                            }
                            if(j<i){
                                S_bond[j]=S_bond[j]+B_matrix_local[j][i];
                            }
                        }
                    }

                    //ノンゼロ成分の行(row1)、列(col1)を記録
                    for(int i=0;i<=H;i=i+1){ //C~H列の探索(X2はリストに含めない)
                        for(int j=i; j<=X; j=j+1){ //C~X行の探索
                            if(S_bond[i]>2 || S_bond[j]>2){ //S原子の結合手が3本以上のときは、反応候補から除外する
                                continue;
                            }
                            if(j==H || j==X){ 
                                for (int k=0 ; k<B_matrix_local[i][j] ; k++){ 
                                    B_nz_local[Total_nz].row1=i;
                                    B_nz_local[Total_nz].col1=j;
                                    Total_nz=Total_nz+1;
                                }
                            }else if (i<C && j>=N){ 
                                if(B_matrix_local[i][j]==1 || B_matrix_local[i][j]==2){//CO分子はカウントしないように
                                    B_nz_local[Total_nz].row1=i;
                                    B_nz_local[Total_nz].col1=j;
                                    Total_nz=Total_nz+1;
                                }
                            }else if(B_matrix_local[i][j]!=0){
                                B_nz_local[Total_nz].row1=i;
                                B_nz_local[Total_nz].col1=j;
                                Total_nz=Total_nz+1;
                            }
                        }
                    }
            

                    //CO二重結合、三重結合成分の行(row1)、列(col1)を記録
                    for(int i=0;i<C;i=i+1){
                        for(int j=N; j<O; j=j+1){
                            if (B_matrix_local[i][j]==2){
                                CO_dou[CO2].row1=i;
                                CO_dou[CO2].col1=j;
                                CO2=CO2+1;
                            }else if (B_matrix_local[i][j]==3){
                                CO_tri[CO3].row1=i;
                                CO_tri[CO3].col1=j;
                                CO3=CO3+1;
                            }
                        }
                    }
            

                    //SO二重結合の行(row1)、列(col1)を記録
                    for(int i=N;i<O;i=i+1){
                        for(int k=O; k<S; k=k+1){
                            if (B_matrix_local[i][k]==2){
                                SO2_mol[SO2].row1=i;    //O原子の番号
                                SO2_mol[SO2].col1=k;    //S原子の番号
                                SO2=SO2+1;
                                sys.bridge[i][k]=1;     //SO二重結合は橋とする
                            }
                        }
                    }
                    
                    //Sの結合数が4または6で、結合数が変わらない反応の行列成分を記録
                    for(int i=0; i<H; i++){
                        int count=0;
                        if(S_bond[i]==4){
                            for(int l=0; l<SO2; l++){
                                if (SO2_mol[l].col1==i){   //O=S(-R)-R' 
                                    if(SO2_mol[l].col1!=count){//ノンゼロ成分の多重カウントを防止
                                        for(int j=0; j<=X; j++){
                                            if(j==H || j==X){
                                                for (int k=0 ; k<B_matrix_local[SO2_mol[l].col1][j] ; k++){
                                                    B_nz_local[Total_nz].row1=SO2_mol[l].col1;    //S原子の番号
                                                    B_nz_local[Total_nz].col1=j;                  //HまたはX原子の番号
                                                    Total_nz=Total_nz+1;
                                                    count=i;
                                                }
                                            }else if(SO2_mol[l].col1<j){
                                                if(B_matrix_local[SO2_mol[l].col1][j]!=0){
                                                    B_nz_local[Total_nz].row1=SO2_mol[l].col1;    //S原子の番号
                                                    B_nz_local[Total_nz].col1=j;                  //RまたはR'原子の番号
                                                    Total_nz=Total_nz+1;
                                                    count=i;
                                                }
                                            }else if(j<SO2_mol[l].col1){
                                                if(B_matrix_local[j][SO2_mol[l].col1]==1){
                                                    B_nz_local[Total_nz].row1=j;                  //RまたはR'原子の番号
                                                    B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                    Total_nz=Total_nz+1;
                                                    count=i;
                                                }else if (B_matrix_local[j][SO2_mol[l].col1]==2){
                                                    if(j<N && O<=j){   //jはO原子を除く
                                                        B_nz_local[Total_nz].row1=j;                  //RまたはR'原子の番号
                                                        B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                        Total_nz=Total_nz+1;
                                                        count=i;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                for(int m=l+1; m<SO2; m++){
                                    if(SO2_mol[l].col1==i && SO2_mol[m].col1==i){   //O=S=O'

                                        B_nz_local[Total_nz].row1=SO2_mol[l].row1;    //O原子の番号
                                        B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                        Total_nz=Total_nz+1;

                                        B_nz_local[Total_nz].row1=SO2_mol[m].row1;    //O'原子の番号
                                        B_nz_local[Total_nz].col1=SO2_mol[m].col1;    //S原子の番号
                                        Total_nz=Total_nz+1;
                                        //break;
                                    }
                                }
                            }
                        }else if (S_bond[i]==6){
                            for(int l=0; l<SO2; l++){
                                for(int m=l+1; m<SO2; m++){
                                    if (SO2_mol[l].col1==i && SO2_mol[m].col1==i){   //O=S(-R)(=O')-R'
                                        if(SO2_mol[l].col1!=count){//ノンゼロ成分の多重カウントを防止
                                            for(int j=0; j<=X; j++){
                                                if(j==H || j==X){
                                                    for (int k=0 ; k<B_matrix_local[SO2_mol[l].col1][j] ; k++){
                                                        B_nz_local[Total_nz].row1=SO2_mol[l].col1;    //S原子の番号
                                                        B_nz_local[Total_nz].col1=j;                  //HまたはX原子の番号
                                                        Total_nz=Total_nz+1;
                                                        count=i;
                                                    }
                                                }else if(SO2_mol[l].col1<j){
                                                    if(B_matrix_local[SO2_mol[l].col1][j]!=0){
                                                        if(S_bond[j]==6){   //Sの結合数が6のときのみ反応候補に加える
                                                            B_nz_local[Total_nz].row1=SO2_mol[l].col1;    //S原子の番号
                                                            B_nz_local[Total_nz].col1=j;                  //RまたはR'原子の番号
                                                            Total_nz=Total_nz+1;
                                                        }
                                                        count=i;
                                                    }
                                                }else if(j<SO2_mol[l].col1 && j<O){
                                                    if(B_matrix_local[j][SO2_mol[l].col1]==1){
                                                        B_nz_local[Total_nz].row1=j;                  //RまたはR'原子の番号
                                                        B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                        Total_nz=Total_nz+1;
                                                        count=i;
                                                    }else if (B_matrix_local[j][SO2_mol[l].col1]==2){
                                                        if(j<N && O<=j){   //jはO原子を除く
                                                            B_nz_local[Total_nz].row1=j;                  //RまたはR'原子の番号
                                                            B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                                            Total_nz=Total_nz+1;
                                                            count=i;
                                                        }
                                                    }
                                                }else{
                                                    count=i;
                                                }
                                            }
                                        }
                                    }
                                    for(int n=m+1; n<SO2; n++){
                                        if(SO2_mol[l].col1==i && SO2_mol[m].col1==i && SO2_mol[n].col1==i){   //O=S(=O')=O"

                                            B_nz_local[Total_nz].row1=SO2_mol[l].row1;    //O原子の番号
                                            B_nz_local[Total_nz].col1=SO2_mol[l].col1;    //S原子の番号
                                            Total_nz=Total_nz+1;

                                            B_nz_local[Total_nz].row1=SO2_mol[m].row1;    //O'原子の番号
                                            B_nz_local[Total_nz].col1=SO2_mol[m].col1;    //S原子の番号
                                            Total_nz=Total_nz+1;

                                            B_nz_local[Total_nz].row1=SO2_mol[n].row1;    //O"原子の番号
                                            B_nz_local[Total_nz].col1=SO2_mol[n].col1;    //S原子の番号
                                            Total_nz=Total_nz+1;
                                            //break;
                                        }
                                    }
                                }
                            }
                        }
                    }//*/

                    //R-O-R'構造の行、列を記録
                    for(int i=N; i<O; i=i+1){
                        for(int j=0; j<=H; j=j+1){
                            for(int k=j+1; k<=X; k=k+1){
                                if(i<j){
                                    if(B_matrix_local[i][j]==1 && B_matrix_local[i][k]==1){
                                        ROR_mol[ROR].row1=i;    //O原子の番号
                                        ROR_mol[ROR].col1=j;    //Rの番号
                                        ROR_mol[ROR].row2=i;    //O原子の番号
                                        ROR_mol[ROR].col2=k;    //R'の番号
                                        ROR=ROR+1;
                                    }
                                }
                                if(j<i && i<k){
                                    if(B_matrix_local[j][i]==1 && B_matrix_local[i][k]==1){
                                        ROR_mol[ROR].row1=i;    //O原子の番号
                                        ROR_mol[ROR].col1=j;    //Rの番号
                                        ROR_mol[ROR].row2=i;    //O原子の番号
                                        ROR_mol[ROR].col2=k;    //R'の番号
                                        ROR=ROR+1;
                                    }
                                }
                                if(k<i){
                                    if(B_matrix_local[j][i]==1 && B_matrix_local[k][i]==1){
                                        ROR_mol[ROR].row1=i;    //O原子の番号
                                        ROR_mol[ROR].col1=j;    //Rの番号
                                        ROR_mol[ROR].row2=i;    //O原子の番号
                                        ROR_mol[ROR].col2=k;    //R'の番号
                                        ROR=ROR+1;
                                    }
                                }
                            }
                        }

                        if(B_matrix_local[i][H]==2){
                            ROR_mol[ROR].row1=i;    //O原子の番号
                            ROR_mol[ROR].col1=H;    //H原子の番号
                            ROR_mol[ROR].row2=i;    //O原子の番号
                            ROR_mol[ROR].col2=H;    //H原子の番号
                            ROR=ROR+1;
                        }
                        if(B_matrix_local[i][X]==2){
                            ROR_mol[ROR].row1=i;    //O原子の番号
                            ROR_mol[ROR].col1=X;    //X原子の番号
                            ROR_mol[ROR].row2=i;    //O原子の番号
                            ROR_mol[ROR].col2=X;    //X原子の番号
                            ROR=ROR+1;
                        }
                    }//*/


                    /*/変化無 NO CHANGE
                    R_nz_local[bb].row1=H, R_nz_local[bb].col1=H; 
                    R_nz_local[bb].row2=H, R_nz_local[bb].col2=H;
                    R_nz_local[bb].row3=H, R_nz_local[bb].col3=H;
                    R_nz_local[bb].row4=H, R_nz_local[bb].col4=H;
                    D_energy_local[bb]=Ea;
                    bb=bb+1;//*/

                    //結合行列のノンゼロリストをもとに、切断させる結合・生成させる結合の組み合わせを作成（反応候補リスト作成）
                    for (int l=0 ; l<Total_nz ; l++){
                    for (int m=l+1 ; m<Total_nz ; m++){ //l番目以降のm番目のノンゼロ成分について
                            if(B_nz_local[l].row1!=B_nz_local[m].row1){  //同じ原子が持つ結合２本は選べない
                                if(B_nz_local[l].col1==H || B_nz_local[l].col1==X){  //R-HまたはR-Xの場合B_nz_local[l].col1とB_nz_local[m].col1は同じでもいい

                                    CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                    R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1;
                                    R_nz_local[bb].row2=B_nz_local[m].row1, R_nz_local[bb].col2=B_nz_local[m].col1;

                                    if(B_nz_local[l].row1<=B_nz_local[m].row1){
                                        R_nz_local[bb].row3=B_nz_local[l].row1, R_nz_local[bb].col3=B_nz_local[m].row1;
                                    }else{
                                        R_nz_local[bb].row3=B_nz_local[m].row1, R_nz_local[bb].col3=B_nz_local[l].row1;
                                    }

                                    if(B_nz_local[l].col1<=B_nz_local[m].col1){
                                        R_nz_local[bb].row4=B_nz_local[l].col1, R_nz_local[bb].col4=B_nz_local[m].col1;
                                    }else{
                                        R_nz_local[bb].row4=B_nz_local[m].col1, R_nz_local[bb].col4=B_nz_local[l].col1;
                                    }
                                    bb = bb + 1;

                                    CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                    R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1;
                                    R_nz_local[bb].row2=B_nz_local[m].row1, R_nz_local[bb].col2=B_nz_local[m].col1;
                                    
                                    if(B_nz_local[l].row1<=B_nz_local[m].col1){
                                        R_nz_local[bb].row3=B_nz_local[l].row1, R_nz_local[bb].col3=B_nz_local[m].col1;
                                    }else{
                                        R_nz_local[bb].row3=B_nz_local[m].col1, R_nz_local[bb].col3=B_nz_local[l].row1;
                                    }

                                    if(B_nz_local[m].row1<=B_nz_local[l].col1){
                                        R_nz_local[bb].row4=B_nz_local[m].row1, R_nz_local[bb].col4=B_nz_local[l].col1;
                                    }else{
                                        R_nz_local[bb].row4=B_nz_local[l].col1, R_nz_local[bb].col4=B_nz_local[m].row1;
                                    }
                                    bb = bb + 1;

                                }else if(B_nz_local[l].col1!=B_nz_local[m].row1 && B_nz_local[l].col1!=B_nz_local[m].col1){ //R1-R2 + R3-R4ならR2はR3,R4と異なる原子でなくてはならない
                                    if(B_nz_local[l].col1<B_nz_local[m].col1){

                                        CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                        R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1;
                                        R_nz_local[bb].row2=B_nz_local[m].row1, R_nz_local[bb].col2=B_nz_local[m].col1;
                                        
                                        if(B_nz_local[l].row1<=B_nz_local[m].row1){
                                            R_nz_local[bb].row3=B_nz_local[l].row1, R_nz_local[bb].col3=B_nz_local[m].row1;
                                        }else{
                                            R_nz_local[bb].row3=B_nz_local[m].row1, R_nz_local[bb].col3=B_nz_local[l].row1;
                                        }

                                        R_nz_local[bb].row4=B_nz_local[l].col1, R_nz_local[bb].col4=B_nz_local[m].col1;

                                        bb = bb + 1;

                                    }else{

                                        CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                        R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1;
                                        R_nz_local[bb].row2=B_nz_local[m].row1, R_nz_local[bb].col2=B_nz_local[m].col1;
                                        
                                        if(B_nz_local[l].row1<=B_nz_local[m].row1){
                                            R_nz_local[bb].row3=B_nz_local[l].row1, R_nz_local[bb].col3=B_nz_local[m].row1;
                                        }else{
                                            R_nz_local[bb].row3=B_nz_local[m].row1, R_nz_local[bb].col3=B_nz_local[l].row1;
                                        }

                                        R_nz_local[bb].row4=B_nz_local[m].col1, R_nz_local[bb].col4=B_nz_local[l].col1;

                                        bb = bb + 1;
                            
                                    }

                                    if(B_nz_local[l].col1<B_nz_local[m].row1){

                                        CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                        R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1;
                                        R_nz_local[bb].row2=B_nz_local[m].row1, R_nz_local[bb].col2=B_nz_local[m].col1;
                                        
                                        if(B_nz_local[l].row1<=B_nz_local[m].col1){
                                            R_nz_local[bb].row3=B_nz_local[l].row1, R_nz_local[bb].col3=B_nz_local[m].col1;
                                        }else{
                                            R_nz_local[bb].row3=B_nz_local[m].col1, R_nz_local[bb].col3=B_nz_local[l].row1;
                                        }

                                        R_nz_local[bb].row4=B_nz_local[l].col1, R_nz_local[bb].col4=B_nz_local[m].row1;

                                        bb = bb + 1;

                                    }else{

                                        CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                        R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1;
                                        R_nz_local[bb].row2=B_nz_local[m].row1, R_nz_local[bb].col2=B_nz_local[m].col1;
                                        
                                        if(B_nz_local[l].row1<=B_nz_local[m].col1){
                                            R_nz_local[bb].row3=B_nz_local[l].row1, R_nz_local[bb].col3=B_nz_local[m].col1;
                                        }else{
                                            R_nz_local[bb].row3=B_nz_local[m].col1, R_nz_local[bb].col3=B_nz_local[l].row1;
                                        }

                                        R_nz_local[bb].row4=B_nz_local[m].row1, R_nz_local[bb].col4=B_nz_local[l].col1;

                                        bb = bb + 1;
                                    }
                        
                                }
                            }
                        }
                    }

                    //CO生成反応 R-CO-R' → R-R' + C三O 
                    for (int i=0 ; i<CO2 ; i++){           
                        for (int k=0 ; k<H ; k++){
                            if (k<CO_dou[i].row1){
                                if (B_matrix_local[k][CO_dou[i].row1]==2){
                                    break;
                                }else if (B_matrix_local[k][CO_dou[i].row1]==1){
                                    for (int l=k+1 ; l<H ; l++){
                                        if (CO_dou[i].row1<l && B_matrix_local[CO_dou[i].row1][l]==1){
                                            CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                            R_nz_local[bb].row1=k             , R_nz_local[bb].col1=CO_dou[i].row1;
                                            R_nz_local[bb].row2=CO_dou[i].row1, R_nz_local[bb].col2=l;
                                            R_nz_local[bb].row3=k             , R_nz_local[bb].col3=l;
                                            R_nz_local[bb].row4=CO_dou[i].row1, R_nz_local[bb].col4=CO_dou[i].col1;
                                            bb=bb+1;
                                            break;
                                        }else if (l<CO_dou[i].row1 && B_matrix_local[l][CO_dou[i].row1]==1){
                                            R_nz_local[bb].row1=k             , R_nz_local[bb].col1=CO_dou[i].row1;
                                            R_nz_local[bb].row2=l             , R_nz_local[bb].col2=CO_dou[i].row1;
                                            R_nz_local[bb].row3=k             , R_nz_local[bb].col3=l;
                                            R_nz_local[bb].row4=CO_dou[i].row1, R_nz_local[bb].col4=CO_dou[i].col1;
                                            bb=bb+1;
                                            break;
                                        }
                                    }
                                }
                            }else if (CO_dou[i].row1<k){
                                if (B_matrix_local[CO_dou[i].row1][k]==2){
                                    break;
                                }else if (B_matrix_local[CO_dou[i].row1][k]==1){
                                    for (int l=k+1 ; l<H ; l++){
                                        if (B_matrix_local[CO_dou[i].row1][l]==1){
                                            CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                            R_nz_local[bb].row1=CO_dou[i].row1, R_nz_local[bb].col1=k;
                                            R_nz_local[bb].row2=CO_dou[i].row1, R_nz_local[bb].col2=l;
                                            R_nz_local[bb].row3=k             , R_nz_local[bb].col3=l;
                                            R_nz_local[bb].row4=CO_dou[i].row1, R_nz_local[bb].col4=CO_dou[i].col1;
                                            bb=bb+1;
                                            break;
                                        }
                                    }
                                }
                            }
                            for(int k=H ; k<=X ; k++){
                                if (B_matrix_local[CO_dou[i].row1][k]==2){
                                    CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                    R_nz_local[bb].row1=CO_dou[i].row1, R_nz_local[bb].col1=k;
                                    R_nz_local[bb].row2=CO_dou[i].row1, R_nz_local[bb].col2=k;
                                    R_nz_local[bb].row3=k,              R_nz_local[bb].col3=k;
                                    R_nz_local[bb].row4=CO_dou[i].row1, R_nz_local[bb].col4=CO_dou[i].col1;
                                    bb=bb+1;
                                    break;
                                }
                            }
                        }
                    }

                    //CO分解反応 CO + R-R' -> R-CO-R'
                    for(int i=0 ; i<CO3 ; i++){ 
                        for (int l=0 ; l<Total_nz ; l++){
                            CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                            R_nz_local[bb].row1=B_nz_local[l].row1,   R_nz_local[bb].col1=B_nz_local[l].col1;
                            R_nz_local[bb].row2=CO_tri[i].row1, R_nz_local[bb].col2=CO_tri[i].col1;

                            if(CO_tri[i].row1<B_nz_local[l].row1){

                                R_nz_local[bb].row3=CO_tri[i].row1, R_nz_local[bb].col3=B_nz_local[l].row1;
                                R_nz_local[bb].row4=CO_tri[i].row1, R_nz_local[bb].col4=B_nz_local[l].col1;
                                bb = bb + 1;

                            }else if (B_nz_local[l].row1<CO_tri[i].row1 && CO_tri[i].row1<B_nz_local[l].col1){

                                R_nz_local[bb].row3=B_nz_local[l].row1,   R_nz_local[bb].col3=CO_tri[i].row1;
                                R_nz_local[bb].row4=CO_tri[i].row1, R_nz_local[bb].col4=B_nz_local[l].col1;
                                bb = bb + 1;

                            }else if (B_nz_local[l].col1<CO_tri[i].row1){

                                R_nz_local[bb].row3=B_nz_local[l].row1,   R_nz_local[bb].col3=CO_tri[i].row1;
                                R_nz_local[bb].row4=B_nz_local[l].col1,   R_nz_local[bb].col4=CO_tri[i].row1;
                                bb = bb + 1;
                            }
                        }
                        for(int i=0 ; i<CO3 ; i++){ //ラジカルとの反応
                            for (int j=0 ; j<X ; j++){
                                for (int k=0 ; k<B_matrix_local[j][X] ; k++){
                                    R_nz_local[bb].row1=j             , R_nz_local[bb].col1=X;
                                        CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                    R_nz_local[bb].row2=CO_tri[i].row1, R_nz_local[bb].col2=CO_tri[i].col1;
                                    R_nz_local[bb].row4=CO_tri[i].row1, R_nz_local[bb].col4=X;
                                    if(CO_tri[i].row1<j){
                                        R_nz_local[bb].row3=CO_tri[i].row1, R_nz_local[bb].col3=j;
                                        bb = bb + 1;
                                    }else{
                                        R_nz_local[bb].row3=j, R_nz_local[bb].col3=CO_tri[i].row1;
                                        bb = bb + 1;
                                    }
                                }
                            }
                        }
                    }

                    int so_1=bb;

                    //Sの結合数が減る反応
                    for(int m=0; m<=X; m++){
                        if(S_bond[m]==4 || S_bond[m]==6){ //Sの結合数が4または6のとき
                            for(int i=0; i<SO2; i++){
                                if(SO2_mol[i].col1==m){ 

                                    //SO二重結合分解反応 S=O + R-R' -> R-O-R' + S
                                    for (int l=0 ; l<Total_nz ; l++){
                                        CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                        if(SO2_mol[i].col1!=B_nz_local[l].row1 && SO2_mol[i].col1!=B_nz_local[l].col1){   //SO二重結合をもつS原子が持つほかの結合は選べない
                                            R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1;
                                            R_nz_local[bb].row2=SO2_mol[i].row1,    R_nz_local[bb].col2=SO2_mol[i].col1; //切断されるO=S結合

                                            if(SO2_mol[i].row1<B_nz_local[l].row1){

                                                R_nz_local[bb].row3=SO2_mol[i].row1,     R_nz_local[bb].col3=B_nz_local[l].row1;
                                                R_nz_local[bb].row4=SO2_mol[i].row1,     R_nz_local[bb].col4=B_nz_local[l].col1;
                                                bb = bb + 1;

                                            }else if (B_nz_local[l].row1<SO2_mol[i].row1 && SO2_mol[i].row1<B_nz_local[l].col1){

                                                R_nz_local[bb].row3=B_nz_local[l].row1,  R_nz_local[bb].col3=SO2_mol[i].row1;
                                                R_nz_local[bb].row4=SO2_mol[i].row1,     R_nz_local[bb].col4=B_nz_local[l].col1;
                                                bb = bb + 1;

                                            }else if (B_nz_local[l].col1<SO2_mol[i].row1){

                                                R_nz_local[bb].row3=B_nz_local[l].row1,  R_nz_local[bb].col3=SO2_mol[i].row1;
                                                R_nz_local[bb].row4=B_nz_local[l].col1,  R_nz_local[bb].col4=SO2_mol[i].row1;
                                                bb = bb + 1;
                                            }
                                        }

                                        //R-S-R'-> R-R' + S
                                        for(int n=l+1; n<Total_nz; n++){
                                            if(B_matrix_local[B_nz_local[l].row1][B_nz_local[l].col1]==1 && B_matrix_local[B_nz_local[n].row1][B_nz_local[n].col1]==1){ //それぞれが単結合であることを確認
                                                if(SO2_mol[i].col1==B_nz_local[l].row1 && SO2_mol[i].col1==B_nz_local[n].row1){
                                                    R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1; //切断されるS-R結合
                                                    R_nz_local[bb].row2=B_nz_local[n].row1, R_nz_local[bb].col2=B_nz_local[n].col1; //切断されるS-R'結合
                                                    R_nz_local[bb].row3=B_nz_local[n].row1, R_nz_local[bb].col3=B_nz_local[n].col1; //2で減りすぎた結合数を補う
                                                    if(B_nz_local[l].col1<=B_nz_local[n].col1){
                                                        R_nz_local[bb].row4=B_nz_local[l].col1, R_nz_local[bb].col4=B_nz_local[n].col1; //生成されるR-R'結合
                                                        bb = bb + 1;
                                                    }else{
                                                        R_nz_local[bb].row4=B_nz_local[n].col1, R_nz_local[bb].col4=B_nz_local[l].col1; //生成されるR'-R結合
                                                        bb = bb + 1;
                                                    }
                                                }else if (SO2_mol[i].col1==B_nz_local[l].row1 && SO2_mol[i].col1==B_nz_local[n].col1){
                                                    R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1; //切断されるS-R結合
                                                    R_nz_local[bb].row2=B_nz_local[n].row1, R_nz_local[bb].col2=B_nz_local[n].col1; //切断されるR'-S結合
                                                    R_nz_local[bb].row3=B_nz_local[n].row1, R_nz_local[bb].col3=B_nz_local[n].col1; //2で減りすぎた結合数を補う
                                                    R_nz_local[bb].row4=B_nz_local[n].row1, R_nz_local[bb].col4=B_nz_local[l].col1; //生成されるR'-R結合
                                                    bb = bb + 1;
                                                }else if (SO2_mol[i].col1==B_nz_local[l].col1 && SO2_mol[i].col1==B_nz_local[n].row1){
                                                    R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1; //切断されるR-S結合
                                                    R_nz_local[bb].row2=B_nz_local[n].row1, R_nz_local[bb].col2=B_nz_local[n].col1; //切断されるS-R'結合
                                                    R_nz_local[bb].row3=B_nz_local[n].row1, R_nz_local[bb].col3=B_nz_local[n].col1; //2で減りすぎた結合数を補う
                                                    R_nz_local[bb].row4=B_nz_local[l].row1, R_nz_local[bb].col4=B_nz_local[n].col1; //生成されるR-R'結合
                                                    bb = bb + 1;
                                                }else if (SO2_mol[i].col1==B_nz_local[l].col1 && SO2_mol[i].col1==B_nz_local[n].col1){
                                                    R_nz_local[bb].row1=B_nz_local[l].row1, R_nz_local[bb].col1=B_nz_local[l].col1; //切断されるR-S結合
                                                    R_nz_local[bb].row2=B_nz_local[n].row1, R_nz_local[bb].col2=B_nz_local[n].col1; //切断されるR'-S結合
                                                    R_nz_local[bb].row3=B_nz_local[n].row1, R_nz_local[bb].col3=B_nz_local[n].col1; //2で減りすぎた結合数を補う
                                                    if(B_nz_local[l].row1<=B_nz_local[n].row1){
                                                        R_nz_local[bb].row4=B_nz_local[l].row1, R_nz_local[bb].col4=B_nz_local[n].row1; //生成されるR-R'結合
                                                        bb = bb + 1;
                                                    }else{
                                                        R_nz_local[bb].row4=B_nz_local[n].row1, R_nz_local[bb].col4=B_nz_local[l].row1; //生成されるR'-R結合
                                                        bb = bb + 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    for(int i=0 ; i<SO2 ; i++){ //ラジカルとの反応
                                        if(SO2_mol[i].col1==m){
                                            for (int j=0 ; j<X ; j++){
                                                for (int k=0 ; k<B_matrix_local[j][X] ; k++){
                                                    CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                                    if(SO2_mol[i].col1!=j){   //SO二重結合をもつS原子が持つほかの結合は選べない
                                                        R_nz_local[bb].row1=j             ,  R_nz_local[bb].col1=X;
                                                        R_nz_local[bb].row2=SO2_mol[i].row1, R_nz_local[bb].col2=SO2_mol[i].col1;
                                                        R_nz_local[bb].row4=SO2_mol[i].row1, R_nz_local[bb].col4=X;
                                                        if(SO2_mol[i].row1<j){
                                                            R_nz_local[bb].row3=SO2_mol[i].row1, R_nz_local[bb].col3=j;
                                                            bb = bb + 1;
                                                        }else{
                                                            R_nz_local[bb].row3=j, R_nz_local[bb].col3=SO2_mol[i].row1;
                                                            bb = bb + 1;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if(B_matrix_local[SO2_mol[i].col1][H]==2){ //S-H結合が2つある場合
                                        R_nz_local[bb].row1=SO2_mol[i].col1, R_nz_local[bb].col1=H;   //3と打ち消しあう
                                        R_nz_local[bb].row2=SO2_mol[i].col1, R_nz_local[bb].col2=H;   //S-H結合の切断
                                        R_nz_local[bb].row3=SO2_mol[i].col1, R_nz_local[bb].col3=H;   //1と打ち消しあう
                                        R_nz_local[bb].row4=H              , R_nz_local[bb].col4=H;   //生成されるH-H結合
                                        bb=bb+1;
                                    }
                                    if(B_matrix_local[SO2_mol[i].col1][X]==2){ //S-X結合が2つある場合
                                        R_nz_local[bb].row1=SO2_mol[i].col1, R_nz_local[bb].col1=X;   //3と打ち消しあう
                                        R_nz_local[bb].row2=SO2_mol[i].col1, R_nz_local[bb].col2=X;   //S-X結合の切断
                                        R_nz_local[bb].row3=SO2_mol[i].col1, R_nz_local[bb].col3=X;   //1と打ち消しあう
                                        R_nz_local[bb].row4=X              , R_nz_local[bb].col4=X;   //生成されるX-X結合
                                        bb=bb+1;
                                    }
                                }
                            }
                        }
                    }

                    int so_2=bb;

                    //Sの結合数が増える反応
                    for(int i=0; i<=X; i++){
                        if(S_bond[i]==2 || S_bond[i]==4){ //Sの結合数が2または4のとき
                            for(int j=0; j<SO2; j++){
                                if(SO2_mol[j].col1==i){

                                    //SO二重結合生成反応 S + R-O-R' -> R-R' + SO
                                    for(int k=0; k<ROR; k++){
                                        CHECK_CAPACITY(bb, R_NZ_CAP, "R_nz_local");
                                        if(SO2_mol[j].col1!=ROR_mol[k].col1 && SO2_mol[j].col1!=ROR_mol[k].col2){   //SO二重結合をもつS原子が持つほかの結合は選べない
                                            if(ROR_mol[k].row1 < ROR_mol[k].col1){
                                                R_nz_local[bb].row1=ROR_mol[k].row1, R_nz_local[bb].col1=ROR_mol[k].col1;   //切断されるO-R結合
                                                R_nz_local[bb].row2=ROR_mol[k].row1, R_nz_local[bb].col2=ROR_mol[k].col2;   //切断されるO-R'結合
                                                R_nz_local[bb].row3=ROR_mol[k].col1, R_nz_local[bb].col3=ROR_mol[k].col2;   //生成されるR-R'結合
                                                R_nz_local[bb].row4=ROR_mol[k].row1, R_nz_local[bb].col4=SO2_mol[j].col1;   //生成されるO=S結合
                                                bb=bb+1;
                                            }
                                            if(ROR_mol[k].col1 < ROR_mol[k].row1 && ROR_mol[k].row1 < ROR_mol[k].col2){
                                                R_nz_local[bb].row1=ROR_mol[k].col1, R_nz_local[bb].col1=ROR_mol[k].row1;   //切断されるR-O結合
                                                R_nz_local[bb].row2=ROR_mol[k].row1, R_nz_local[bb].col2=ROR_mol[k].col2;   //切断されるO-R'結合
                                                R_nz_local[bb].row3=ROR_mol[k].col1, R_nz_local[bb].col3=ROR_mol[k].col2;   //生成されるR-R'結合
                                                R_nz_local[bb].row4=ROR_mol[k].row1, R_nz_local[bb].col4=SO2_mol[j].col1;   //生成されるO=S結合
                                                bb=bb+1;
                                            }
                                            if(ROR_mol[k].row1 > ROR_mol[k].col2){
                                                R_nz_local[bb].row1=ROR_mol[k].col1, R_nz_local[bb].col1=ROR_mol[k].row1;   //切断されるR-O結合
                                                R_nz_local[bb].row2=ROR_mol[k].col2, R_nz_local[bb].col2=ROR_mol[k].row1;   //切断されるR'-O結合
                                                R_nz_local[bb].row3=ROR_mol[k].col1, R_nz_local[bb].col3=ROR_mol[k].col2;   //生成されるR-R'結合
                                                R_nz_local[bb].row4=ROR_mol[k].row1, R_nz_local[bb].col4=SO2_mol[j].col1;   //生成されるO=S結合
                                                bb=bb+1;
                                            }
                                        }

                                        //S + R-R' -> R-S-R'
                                        int count_2=0; //多重カウントを防止するためのフラグ
                                        int count_4=0; //多重カウントを防止するためのフラグ
                                        for(int l=j+1; l<SO2; l++){
                                            if(S_bond[j]==2 && count_2!=l){   //多重カウントを防止
                                                for(int n=0; n<Total_nz; n++){
                                                    if(SO2_mol[j].col1!=B_nz_local[n].row1 && SO2_mol[j].col1!=B_nz_local[n].col1){   //SO二重結合をもつS原子が持つほかの結合は選べない
                                                        R_nz_local[bb].row1=B_nz_local[n].row1, R_nz_local[bb].col1=B_nz_local[n].col1; //切断されるR-R'結合
                                                        if(SO2_mol[j].col1<B_nz_local[n].row1){
                                                            R_nz_local[bb].row2=SO2_mol[j].col1,    R_nz_local[bb].col2=B_nz_local[n].col1; //4で増えすぎた結合数を補う
                                                            R_nz_local[bb].row3=SO2_mol[j].col1,    R_nz_local[bb].col3=B_nz_local[n].row1; //生成されるS-R結合
                                                            R_nz_local[bb].row4=SO2_mol[j].col1,    R_nz_local[bb].col4=B_nz_local[n].col1; //生成されるS-R'結合
                                                            bb = bb + 1;
                                                        }else if (B_nz_local[n].row1<SO2_mol[j].col1 && SO2_mol[j].col1<B_nz_local[n].col1){
                                                            R_nz_local[bb].row2=SO2_mol[j].col1,    R_nz_local[bb].col2=B_nz_local[n].col1;    //4で増えすぎた結合数を補う
                                                            R_nz_local[bb].row3=B_nz_local[n].row1, R_nz_local[bb].col3=SO2_mol[j].col1; //生成されるR-S結合
                                                            R_nz_local[bb].row4=SO2_mol[j].col1,    R_nz_local[bb].col4=B_nz_local[n].col1;    //生成されるS-R'結合
                                                            bb = bb + 1;
                                                        }else if (B_nz_local[n].col1<SO2_mol[j].col1){
                                                            R_nz_local[bb].row2=B_nz_local[n].col1, R_nz_local[bb].col2=SO2_mol[j].col1; //4で増えすぎた結合数を補う
                                                            R_nz_local[bb].row3=B_nz_local[n].row1, R_nz_local[bb].col3=SO2_mol[j].col1; //生成されるR-S結合
                                                            R_nz_local[bb].row4=B_nz_local[n].col1, R_nz_local[bb].col4=SO2_mol[j].col1; //生成されるR'-S結合
                                                            bb = bb + 1;
                                                        }
                                                    }
                                                    count_2 = l;
                                                }
                                            }else if (S_bond[j]==4 && count_4!=l){   //多重カウントを防止
                                                if(SO2_mol[j].col1==i && SO2_mol[l].col1==i){   //S原子はSO二重結合以外をもってはいけない
                                                    for(int n=0; n<Total_nz; n++){
                                                        if(SO2_mol[j].col1!=B_nz_local[n].row1 && SO2_mol[j].col1!=B_nz_local[n].col1){   //SO二重結合をもつS原子が持つほかの結合は選べない
                                                            R_nz_local[bb].row1=B_nz_local[n].row1, R_nz_local[bb].col1=B_nz_local[n].col1; //切断されるR-R'結合
                                                            if(SO2_mol[j].col1<B_nz_local[n].row1){
                                                                R_nz_local[bb].row2=SO2_mol[j].col1,    R_nz_local[bb].col2=B_nz_local[n].col1; //4で増えすぎた結合数を補う
                                                                R_nz_local[bb].row3=SO2_mol[j].col1,    R_nz_local[bb].col3=B_nz_local[n].row1; //生成されるS-R結合
                                                                R_nz_local[bb].row4=SO2_mol[j].col1,    R_nz_local[bb].col4=B_nz_local[n].col1; //生成されるS-R'結合
                                                                bb = bb + 1;
                                                            }else if (B_nz_local[n].row1<SO2_mol[j].col1 && SO2_mol[j].col1<B_nz_local[n].col1){
                                                                R_nz_local[bb].row2=SO2_mol[j].col1,    R_nz_local[bb].col2=B_nz_local[n].col1;    //4で増えすぎた結合数を補う
                                                                R_nz_local[bb].row3=B_nz_local[n].row1, R_nz_local[bb].col3=SO2_mol[j].col1; //生成されるR-S結合
                                                                R_nz_local[bb].row4=SO2_mol[j].col1,    R_nz_local[bb].col4=B_nz_local[n].col1;    //生成されるS-R'結合
                                                                bb = bb + 1;
                                                            }else if (B_nz_local[n].col1<SO2_mol[j].col1){
                                                                R_nz_local[bb].row2=B_nz_local[n].col1, R_nz_local[bb].col2=SO2_mol[j].col1; //4で増えすぎた結合数を補う
                                                                R_nz_local[bb].row3=B_nz_local[n].row1, R_nz_local[bb].col3=SO2_mol[j].col1; //生成されるR-S結合
                                                                R_nz_local[bb].row4=B_nz_local[n].col1, R_nz_local[bb].col4=SO2_mol[j].col1; //生成されるR'-S結合
                                                                bb = bb + 1;
                                                            }
                                                        }
                                                        count_4 = l;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    int so_3=bb;

                    /*光解離反応 R-R' + X2 -> R-X + R'-X
                    for (int l=0 ; l<Total_nz ; l++){
                        for (int j=0 ; j<B_matrix[X][X] ; j++){
                            if(B_nz[l].col1!=X){
                                R_nz_local[bb].row1=B_nz[l].row1, R_nz_local[bb].col1=B_nz[l].col1;
                                R_nz_local[bb].row2=X,            R_nz_local[bb].col2=X;
                                R_nz_local[bb].row3=B_nz[l].row1, R_nz_local[bb].col3=X;
                                R_nz_local[bb].row4=B_nz[l].col1, R_nz_local[bb].col4=X;
                                D_energy[bb]=0;//光解離反応は活性化エネルギーなし
                                bb=bb+1;
                            }
                        }
                    }*/

                    for (int i=0 ; i<bb ; i++){
                        dH_local[i]=B_energy[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]][R_nz_local[i].row1][R_nz_local[i].col1]//切る
                                    -B_energy[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]-1][R_nz_local[i].row1][R_nz_local[i].col1]

                                    +B_energy[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]][R_nz_local[i].row2][R_nz_local[i].col2]//切る
                                    -B_energy[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]-1][R_nz_local[i].row2][R_nz_local[i].col2]

                                    +B_energy[B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3]][R_nz_local[i].row3][R_nz_local[i].col3]
                                    -B_energy[B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3]+1][R_nz_local[i].row3][R_nz_local[i].col3]//結合

                                    +B_energy[B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4]][R_nz_local[i].row4][R_nz_local[i].col4]
                                    -B_energy[B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4]+1][R_nz_local[i].row4][R_nz_local[i].col4];//結合


                        /*/並進運動によるエントロピー（気体の場合）
                        if(atom[R_nz_local[i].row1].molecule != atom[R_nz_local[i].row2].molecule){//2つの結合が異なる分子のとき　二次反応の場合、反応速度にRTNをかける必要があるため、ここで考慮する

                            if(sys.bridge[R_nz_local[i].row1][R_nz_local[i].col1]==1 && sys.bridge[R_nz_local[i].row2][R_nz_local[i].col2]==1){//2つの結合が橋結合のとき
                                dS_local[i]=R*1.5*log(pow((molecule[atom[R_nz_local[i].row1].molecule].mw + molecule[atom[R_nz_local[i].row2].molecule].mw)/2.0, 2)
                                            / (molecule[atom[R_nz_local[i].row1].molecule].mw * molecule[atom[R_nz_local[i].row2].molecule].mw)) -(log(R*T/current_molecule))/(alpha*T);
                            }
                            else{//どちらかが橋結合でないとき
                                dS_local[i]=R*1.5*log((molecule[atom[R_nz_local[i].row1].molecule].mw + molecule[atom[R_nz_local[i].row2].molecule].mw)
                                            / (molecule[atom[R_nz_local[i].row1].molecule].mw * molecule[atom[R_nz_local[i].row2].molecule].mw))-(2.5*log(T/300)+13) -(log(R*T/current_molecule))/(alpha*T);
                            }
                        }
                        else{//2つの結合が同じ分子のとき

                            if(sys.bridge[R_nz_local[i].row1][R_nz_local[i].col1]==0 || sys.bridge[R_nz_local[i].row2][R_nz_local[i].col2]==0){//どちらかが橋結合でないとき
                                dS_local[i]=0.;
                            }
                            else{//2つの結合が橋結合のとき

                                for(int i=0; i<H; i++){
                                    test_sys.visited[i]=0; //分子の探索済みフラグの初期化
                                }

                                B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] + 1;//結合したと仮定
                                B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] + 1;//結合したと仮定

                                DFS(&test_sys, test_atom, test_mol, R_nz_local[i].row3, R_nz_local[i].row3, B_matrix_local, H);
                                DFS(&test_sys, test_atom, test_mol, R_nz_local[i].row4, R_nz_local[i].row4, B_matrix_local, H);

                                if(test_sys.bridge[R_nz_local[i].row3][R_nz_local[i].col3]==1 && test_sys.bridge[R_nz_local[i].row4][R_nz_local[i].col4]==1){//生成された2つの結合が橋結合のとき
                                    dS_local[i]=0.;
                                }
                                else{//生成された2つの結合のうち、 どちらかが橋結合でないとき
                                    dS_local[i]=R*1.5*log(pow((molecule[atom[R_nz_local[i].row1].molecule].mw + molecule[atom[R_nz_local[i].row2].molecule].mw)/2.0, 2)
                                                / molecule[atom[R_nz_local[i].row1].molecule].mw)+(R*2.5*log(T/300)+R*13);
                                }

                                B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] - 1;//元に戻す
                                B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] - 1;//元に戻す
                            
                            }
                        }//*/

                        //振動によるエントロピー（固体の場合）
                        if(atom[R_nz_local[i].row1].molecule != atom[R_nz_local[i].row2].molecule){//2つの結合が異なる分子のとき　二次反応の場合、反応速度にRTNをかける必要があるため、ここで考慮する

                            if(sys.bridge[R_nz_local[i].row1][R_nz_local[i].col1]==1 && sys.bridge[R_nz_local[i].row2][R_nz_local[i].col2]==1){//2つの結合が橋結合のとき
                                int r1_T = molecule[atom[R_nz_local[i].row1].molecule].theta;   //反応物1の振動特性温度
                                int r1_N = molecule[atom[R_nz_local[i].row2].molecule].nAtom;   //反応物1の原子数
                                int r2_T = molecule[atom[R_nz_local[i].row1].molecule].theta;   //反応物2の振動特性温度
                                int r2_N = molecule[atom[R_nz_local[i].row2].molecule].nAtom;   //反応物2の原子数
                                int p_T = (molecule[atom[R_nz_local[i].row1].molecule].theta + molecule[atom[R_nz_local[i].row2].molecule].theta)/2.0;  //生成物の振動特性温度
                                int p_N = (molecule[atom[R_nz_local[i].row1].molecule].nAtom + molecule[atom[R_nz_local[i].row2].molecule].nAtom)/2.0;  //生成物の原子数

                                dS_local[i]=R*(2*max(3*p_N-6, 1)*((p_T/T)/(exp(p_T/T)-1.0) - log(1.0 - exp(-p_T/T)))
                                                - max(3*p_N-6, 1)*((r1_T/T)/(exp(r1_T/T)-1.0) - log(1.0 - exp(-r1_T/T)))
                                                - max(3*p_N-6, 1)*((r2_T/T)/(exp(r2_T/T)-1.0) - log(1.0 - exp(-r2_T/T)))) -(log(R*T/current_molecule))/(alpha*T);
                            }
                            else{//どちらかが橋結合でないとき
                                int r1_T = molecule[atom[R_nz_local[i].row1].molecule].theta;   //反応物1の振動特性温度
                                int r1_N = molecule[atom[R_nz_local[i].row2].molecule].nAtom;   //反応物1の原子数
                                int r2_T = molecule[atom[R_nz_local[i].row1].molecule].theta;   //反応物2の振動特性温度
                                int r2_N = molecule[atom[R_nz_local[i].row2].molecule].nAtom;   //反応物2の原子数
                                int p_T = (molecule[atom[R_nz_local[i].row1].molecule].theta + molecule[atom[R_nz_local[i].row2].molecule].theta)/2.0;  //生成物の振動特性温度
                                int p_N = (molecule[atom[R_nz_local[i].row1].molecule].nAtom + molecule[atom[R_nz_local[i].row2].molecule].nAtom);      //生成物の原子数

                                dS_local[i]=R*(max(3*p_N-6, 1)*((p_T/T)/(exp(p_T/T)-1.0) - log(1.0 - exp(-p_T/T)))
                                                - max(3*p_N-6, 1)*((r1_T/T)/(exp(r1_T/T)-1.0) - log(1.0 - exp(-r1_T/T)))
                                                - max(3*p_N-6, 1)*((r2_T/T)/(exp(r2_T/T)-1.0) - log(1.0 - exp(-r2_T/T)))) -(log(R*T/current_molecule))/(alpha*T);
                            }
                        }
                        else{//2つの結合が同じ分子のとき

                            if(sys.bridge[R_nz_local[i].row1][R_nz_local[i].col1]==0 || sys.bridge[R_nz_local[i].row2][R_nz_local[i].col2]==0){//どちらかが橋結合でないとき
                                dS_local[i]=0.;
                            }
                            else{//2つの結合が橋結合のとき

                                for(int i=0; i<H; i++){
                                    test_sys.visited[i]=0; //分子の探索済みフラグの初期化
                                }

                                B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] + 1;//結合したと仮定
                                B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] + 1;//結合したと仮定

                                DFS(&test_sys, test_atom, test_mol, R_nz_local[i].row3, R_nz_local[i].row3, B_matrix_local, H);
                                DFS(&test_sys, test_atom, test_mol, R_nz_local[i].row4, R_nz_local[i].row4, B_matrix_local, H);

                                if(test_sys.bridge[R_nz_local[i].row3][R_nz_local[i].col3]==1 && test_sys.bridge[R_nz_local[i].row4][R_nz_local[i].col4]==1){//生成された2つの結合が橋結合のとき
                                    dS_local[i]=0.;
                                }
                                else{//生成された2つの結合のうち、 どちらかが橋結合でないとき
                                    int r_T = molecule[atom[R_nz_local[i].row1].molecule].theta;        //反応物2の原子数
                                    int r_N = molecule[atom[R_nz_local[i].row1].molecule].nAtom;        //反応物2の原子数
                                    int p_T = molecule[atom[R_nz_local[i].row1].molecule].theta;        //生成物の振動特性温度
                                    int p_N = molecule[atom[R_nz_local[i].row1].molecule].nAtom/2.0;    //生成物の原子数

                                    dS_local[i]=R*(2*max(3*p_N-6, 1)*((p_T/T)/(exp(p_T/T)-1.0) - log(1.0 - exp(-p_T/T)))
                                                    - max(3*p_N-6, 1)*((r_T/T)/(exp(r_T/T)-1.0) - log(1.0 - exp(-r_T/T))));
                                }

                                B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] - 1;//元に戻す
                                B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] - 1;//元に戻す
                            
                            }
                        }//*/
                        

                        if(i>=so_1 && i<so_2){ //SO2分解反応のとき
                            dH_local[i]=B_energy[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]][R_nz_local[i].row1][R_nz_local[i].col1]//切る
                                    -B_energy[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]-1][R_nz_local[i].row1][R_nz_local[i].col1]

                                    +B_energy[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]][R_nz_local[i].row2][R_nz_local[i].col2]//切る
                                    -B_energy[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]-2][R_nz_local[i].row2][R_nz_local[i].col2]

                                    +B_energy[B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3]][R_nz_local[i].row3][R_nz_local[i].col3]
                                    -B_energy[B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3]+1][R_nz_local[i].row3][R_nz_local[i].col3]//結合

                                    +B_energy[B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4]][R_nz_local[i].row4][R_nz_local[i].col4]
                                    -B_energy[B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4]+1][R_nz_local[i].row4][R_nz_local[i].col4];//結合
                        }

                        if(i>=so_2 && i<so_3){ //SO2生成反応のとき
                            dH_local[i]=B_energy[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]][R_nz_local[i].row1][R_nz_local[i].col1]//切る
                                    -B_energy[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]-1][R_nz_local[i].row1][R_nz_local[i].col1]

                                    +B_energy[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]][R_nz_local[i].row2][R_nz_local[i].col2]//切る
                                    -B_energy[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]-1][R_nz_local[i].row2][R_nz_local[i].col2]

                                    +B_energy[B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3]][R_nz_local[i].row3][R_nz_local[i].col3]
                                    -B_energy[B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3]+1][R_nz_local[i].row3][R_nz_local[i].col3]//結合

                                    +B_energy[B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4]][R_nz_local[i].row4][R_nz_local[i].col4]
                                    -B_energy[B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4]+2][R_nz_local[i].row4][R_nz_local[i].col4];//結合
                        }

                        if(R_nz_local[i].row1==X && R_nz_local[i].col1==X){//光解離反応
                            dH_local[i]=100.;//100は適当な値です
                        }else if(R_nz_local[i].row2==X && R_nz_local[i].col2==X){//光解離反応
                            dH_local[i]=100.;//100は適当な値です
                        }else if(R_nz_local[i].row4==X && R_nz_local[i].col4==X){//ラジカル-ラジカル反応
                            dH_local[i]=dH_local[i]+B_energy[1][X][X];
                            D_energy_local[i]=dH_local[i];
                        }else if(R_nz_local[i].row3==X && R_nz_local[i].col3==X){//ラジカル-ラジカル反応
                            dH_local[i]=dH_local[i]+B_energy[1][X][X];
                            D_energy_local[i]=dH_local[i];
                        }else{
                            D_energy_local[i]=alpha*dH_local[i]-alpha*T*dS_local[i]+beta;
                        }

                        if(D_energy_local[i]<0){//活性化エネルギーはマイナスにならないように
                            D_energy_local[i]=0;
                        }
                    }

                    //HX+HX->H2+X2
                    for(int i=0 ; i<(B_matrix_local[H][X]*(B_matrix_local[H][X]-1))/2 ; i++){
                        R_nz_local[bb].row1=H, R_nz_local[bb].col1=X;
                        R_nz_local[bb].row2=H, R_nz_local[bb].col2=X;
                        R_nz_local[bb].row3=H, R_nz_local[bb].col3=H;
                        R_nz_local[bb].row4=X, R_nz_local[bb].col4=X;
                        D_energy_local[bb]=0;//activation E
                        dH_local[bb]=-B_energy[1][H][H];
                        bb=bb+1;
                    }

                    int Total_R1=bb;

                    //C三O+XOH->HO-CX=O
                    for(int i=0 ; i<CO3 ; i++){ 
                        for (int j=N ; j<O ; j++){
                            if(B_matrix_local[j][H]==1 && B_matrix_local[j][X]==1){
                                R_nz_local[bb].row1=j             , R_nz_local[bb].col1=X;
                                R_nz_local[bb].row2=CO_tri[i].row1, R_nz_local[bb].col2=CO_tri[i].col1;
                                R_nz_local[bb].row3=CO_tri[i].row1, R_nz_local[bb].col3=j;
                                R_nz_local[bb].row4=CO_tri[i].row1, R_nz_local[bb].col4=X;
                                D_energy_local[bb]=100;//activation E
                                dH_local[bb]=-31.;
                                bb = bb + 1;
                            }
                        }
                    }

                    int Total_R=bb;

                    //立体因子Pの計算

                    int n=0;    //結合する原子数

                    for (int i=0 ; i<Total_R ; i++){
                        n=0;
                        if(R_nz_local[i].row1==H && R_nz_local[i].col1==H){//水素解離反応

                            for(int j=0 ; j<H ; j++){
                                if(j!=R_nz_local[i].col2){
                                    if(j<R_nz_local[i].row2 && B_matrix_local[j][R_nz_local[i].row2]>0){
                                        n=n+1;
                                    }else if(j>R_nz_local[i].row2 && B_matrix_local[R_nz_local[i].row2][j]>0){
                                        n=n+1;
                                    }
                                }else if(j!=R_nz_local[i].row2){
                                    if(j<R_nz_local[i].col2 && B_matrix_local[j][R_nz_local[i].col2]>0){
                                        n=n+1;
                                    }else if(j>R_nz_local[i].col2 && B_matrix_local[R_nz_local[i].col2][j]>0){
                                        n=n+1;
                                    }
                                }
                            }
                            n=n+B_matrix_local[R_nz_local[i].row2][H] + B_matrix_local[R_nz_local[i].row2][X] + B_matrix_local[R_nz_local[i].col2][H] + B_matrix_local[R_nz_local[i].col2][X];

                            P_local[i] = (B_length[1][H][H] + B_length[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]][R_nz_local[i].row2][R_nz_local[i].col2])
                                        *(B_length[1][H][H] + B_length[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]][R_nz_local[i].row2][R_nz_local[i].col2])
                                        /**((18-n)/18.0)*/;    //立体構造を表現　　18はまわりをとりまく空間に充填できる原子の数（六方最密構造を仮定）

                        }else if(R_nz_local[i].row2==H && R_nz_local[i].col2==H){//水素解離反応

                            for(int j=0 ; j<H ; j++){
                                if(j!=R_nz_local[i].col1){
                                    if(j<R_nz_local[i].row1 && B_matrix_local[j][R_nz_local[i].row1]>0){
                                        n=n+1;
                                    }else if(j>R_nz_local[i].row1 && B_matrix_local[R_nz_local[i].row1][j]>0){
                                        n=n+1;
                                    }
                                }else if(j!=R_nz_local[i].row1){
                                    if(j<R_nz_local[i].col1 && B_matrix_local[j][R_nz_local[i].col1]>0){
                                        n=n+1;
                                    }else if(j>R_nz_local[i].col1 && B_matrix_local[R_nz_local[i].col1][j]>0){
                                        n=n+1;
                                    }
                                }
                            }
                            n=n+B_matrix_local[R_nz_local[i].row1][H] + B_matrix_local[R_nz_local[i].row1][X] + B_matrix_local[R_nz_local[i].col1][H] + B_matrix_local[R_nz_local[i].col1][X];

                            P_local[i] = (B_length[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]][R_nz_local[i].row1][R_nz_local[i].col1] + B_length[1][H][H])
                                        *(B_length[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]][R_nz_local[i].row1][R_nz_local[i].col1] + B_length[1][H][H])
                                        /**((18-n)/18.0)*/;
                        }else{

                            for(int j=0 ; j<H ; j++){
                                if(j!=R_nz_local[i].col1){
                                    if(j<R_nz_local[i].row1 && B_matrix_local[j][R_nz_local[i].row1]>0){
                                        n=n+1;
                                    }else if(j>R_nz_local[i].row1 && B_matrix_local[R_nz_local[i].row1][j]>0){
                                        n=n+1;
                                    }
                                }else if(j!=R_nz_local[i].row1){
                                    if(j<R_nz_local[i].col1 && B_matrix_local[j][R_nz_local[i].col1]>0){
                                        n=n+1;
                                    }else if(j>R_nz_local[i].col1 && B_matrix_local[R_nz_local[i].col1][j]>0){
                                        n=n+1;
                                    }
                                }
                            }
                            n=n+B_matrix_local[R_nz_local[i].row1][H] + B_matrix_local[R_nz_local[i].row1][X] + B_matrix_local[R_nz_local[i].col1][H] + B_matrix_local[R_nz_local[i].col1][X];

                            for(int j=0 ; j<H ; j++){
                                if(j!=R_nz_local[i].col2){
                                    if(j<R_nz_local[i].row2 && B_matrix_local[j][R_nz_local[i].row2]>0){
                                        n=n+1;
                                    }else if(j>R_nz_local[i].row2 && B_matrix_local[R_nz_local[i].row2][j]>0){
                                        n=n+1;
                                    }
                                }else if(j!=R_nz_local[i].row2){
                                    if(j<R_nz_local[i].col2 && B_matrix_local[j][R_nz_local[i].col2]>0){
                                        n=n+1;
                                    }else if(j>R_nz_local[i].col2 && B_matrix_local[R_nz_local[i].col2][j]>0){
                                        n=n+1;
                                    }
                                }
                            }
                            n=n+B_matrix_local[R_nz_local[i].row2][H] + B_matrix_local[R_nz_local[i].row2][X] + B_matrix_local[R_nz_local[i].col2][H] + B_matrix_local[R_nz_local[i].col2][X];

                            P_local[i] = (B_length[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]][R_nz_local[i].row1][R_nz_local[i].col1]
                                          +B_length[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]][R_nz_local[i].row2][R_nz_local[i].col2])
                                        *(B_length[B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1]][R_nz_local[i].row1][R_nz_local[i].col1]
                                          +B_length[B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2]][R_nz_local[i].row2][R_nz_local[i].col2])
                                        /**((18-n)/18.0)*/;
                        }
                    }//*/

                    /*/
                    for(int i=0 ; i<Total_R ; i++){
                        P_local[i] = 1.0; //立体因子は1.0に固定
                    }//*/

                    //max=立体因子の最大値
                    double max = P_local[0];
                    for (int i=1 ; i<Total_R ; i++){
                        if( P_local[i] > max ){
                            max=P_local[i];
                        }
                    }

                    P_local[0]=max;  //立体因子の最大値をP_local[0]に格納//*/

                    //min=活性化エネルギーの最小値
                    double min = D_energy_local[0];
                    for (int i=1 ; i<Total_R ; i++){
                        if ( D_energy_local[i] < min ){
                            min = D_energy_local[i];
                        }
                        if(min==0.){
                            break;
                        }
                    }

                    //double kTh=6.25e12; //係数
                    double kTh=1.0; //kThの影響を考慮しない

                    //反応速度の和
                    double Total_exp=0;
                    for (int i=0 ; i<Total_R ; i++){
                        if ((D_energy_local[i]-min)<590.*(T/100.)){//underflow防止
                            Total_exp = Total_exp + ((P_local[i]/max) * kTh * exp(-(D_energy_local[i]-min)/(R*T)));
                        }
                    }

                    //r0_local = fmod(ra*r0_local,rm); //乱数

                    //double Prob_reaction=Prob_reaction=1.0-exp(-Total_exp*dt);  //反応が起こる確率

                    //printf("STEP=%d Total_exp=%e prob=%lf bb=%d\n", J, Total_exp, Prob_reaction, Total_R);
                    //if(J==STEP-1) printf("P(%d)=%lf\n", 100, 1.0-exp(-100.0));

                    //if(r0_local/rm < Prob_reaction){//反応が起こる条件
                        //反応が起こるとき

                        //反応の選択

                        //ランダムに反応選択 → E行列（生成物行列）生成
                        int re;
                        double ex ;
                        ex = 0 ;
                        r0_local = fmod(ra*r0_local,rm); //乱数
                        for (int i=0 ; i<Total_R ; i++ ){
                            if ((D_energy_local[i]-min)<590.*(T/100.)){
                                if (ex <= r0_local/rm && r0_local/rm < ex + ((P_local[i]/max) * kTh * exp(-(D_energy_local[i]-min)/(R*T))/Total_exp)){
                                    re=i;//選ばれた反応の番号

                                    if(i>=so_1 && i<so_2){ //SO2分解反応のとき
                                        B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] + 1;//made
                                        B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] + 1;//made
                                        B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] = B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] - 1;//broken
                                        B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] = B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] - 2;//broken
                                        break;
                                    }

                                    if(i>=so_2 && i<so_3){ //SO2生成反応のとき
                                        B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] + 1;//made
                                        B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] + 2;//made
                                        B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] = B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] - 1;//broken
                                        B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] = B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] - 1;//broken
                                        break;
                                    }

                                        B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] + 1;//made
                                        B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] + 1;//made
                                        B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] = B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] - 1;//broken
                                        B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] = B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] - 1;//broken
                                    if(i>=Total_R1){
                                        B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] = B_matrix_local[R_nz_local[i].row3][R_nz_local[i].col3] + 1;//C=Omade
                                        B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] = B_matrix_local[R_nz_local[i].row4][R_nz_local[i].col4] - 1;//C-Xbroken
                                        B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] = B_matrix_local[R_nz_local[i].row1][R_nz_local[i].col1] - 1;//O-Hbroken
                                        B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] = B_matrix_local[R_nz_local[i].row2][R_nz_local[i].col2] + 1;//H-Xmade
                                    }
                                    break;
                                }
                                ex = ex + (P_local[i]/max) * kTh * exp(-(D_energy_local[i]-min)/(R*T))/Total_exp;
                            }
                        }

                        if(re!=0){
                            if(J==PostUV-1){
                                collision_th_UV++; //衝突の回数をカウント
                            }else if((STEP!=PostUV) && (J==STEP-1)){
                                collision_th_PostUV++; //衝突の回数をカウント
                            }else if(J%interval==0){
                                collision_th[J/interval]=collision_th[J/interval]+1; //衝突の回数をカウント
                            }
                        }
                        
                    //}
                }

                //--------------------------------衝突による反応候補の選択終了--------------------------------



                //--------------------------------分子種の判別--------------------------------
                /*/CO, CO2
                for(int i=0;i<C;i=i+1){
                    for(int j=N; j<O; j=j+1){
                        if(B_matrix_local[i][j]==3){
                            CO_abundance[J]=CO_abundance[J]+1;
                            break;
                        }else if(B_matrix_local[i][j]==2){
                            for (int k=j+1; k<O; k++){
                                if (B_matrix_local[i][k]==2){
                                    CO2_abundance[J]=CO2_abundance[J]+1;
                                    break;
                                }
                            }

                        }
                    }
                }
                //H2
                H2_abundance[J]=H2_abundance[J]+B_matrix_local[H][H];

                //Water
                for(int i=N ; i<O ; i++){
                    if(B_matrix_local[i][H]==2){
                        H2O_abundance[J]=H2O_abundance[J]+1;
                    }
                }//*/

                int run=TRIAL;
                int step=STEP;

                //アミノ酸
                if(J%interval==0 || J==PostUV-1 || J==STEP-1){ //計算コスト削減のため、UV phaseでは20回に一回だけ探索
                    //アミノ基(C-NH2)
                    int amino[100];
                    //カルボキシ基
                    int carboxy[100];

                    int AM=0;
                    int CA=0;

                    for (int i=0 ; i<C ; i++){


                        for (int j=C ; j<N ; j++){
                            if(B_matrix_local[i][j]==1 && B_matrix_local[j][H]==2){
                                amino[AM]=i;
                                AM=AM+1;
                                break;
                            }
                        }

                        for (int j=N ; j<O ; j++){
                            if (B_matrix_local[i][j]==1 && B_matrix_local[j][H]==1){//ヒドロキシ基
                                for (int l=N ; l<O ; l++){
                                    if (B_matrix_local[i][l]==2){//カルボキシル基
                                        carboxy[CA]=i;
                                        CA=CA+1;
                                        break;
                                    }
                                }
                            }
                        }
                    }


                    //アミノ酸
                    for (int i=0 ; i<AM ; i++){
                        for (int j=0 ; j<CA ; j++){
                            if (amino[i]<carboxy[j]){
                                if (B_matrix_local[amino[i]][carboxy[j]]==1 && B_matrix_local[amino[i]][H]>=1){
                                    if(J==PostUV-1){
                                        aacid_abundance_UV=aacid_abundance_UV+1;
                                    }else if((STEP!=PostUV) && (J==STEP-1)){
                                        aacid_abundance_PostUV=aacid_abundance_PostUV+1;
                                    }else{
                                        aacid_abundance[J/interval]=aacid_abundance[J/interval]+1;
                                    }
                                }
                            }else if(carboxy[j]<amino[i]){
                                if (B_matrix_local[carboxy[j]][amino[i]]==1 && B_matrix_local[amino[i]][H]>=1){
                                    if(J==PostUV-1){
                                        aacid_abundance_UV=aacid_abundance_UV+1;
                                    }else if((STEP!=PostUV) && (J==STEP-1)){
                                        aacid_abundance_PostUV=aacid_abundance_PostUV+1;
                                    }else{
                                        aacid_abundance[J/interval]=aacid_abundance[J/interval]+1;
                                    }
                                }
                            }
                        }
                    }

                    // Build residues list for peptide detection
                    int residueCount = 0;
                    struct Residue residues[100];
                    int usedPair[80][80];
                    for (int ui=0; ui<80; ui++){
                        for (int uj=0; uj<80; uj++) usedPair[ui][uj]=0;
                    }

                    //amine[], carbonyl[]をアミノ酸残基用にカウント
                    //C-NH- or C-NH2
                    int amine[100];
                    //C=O
                    int carbonyl[100];

                    int am=0;
                    int ca=0;
                    for (int i=0 ; i<C ; i++){

                        for (int j=C ; j<N ; j++){
                            if(B_matrix_local[i][j]==1 && B_matrix_local[j][H]>=1){
                                amine[am]=i;
                                am=am+1;
                                break;
                            }
                        }

                        for (int j=N ; j<O ; j++){                        
                            for (int l=N ; l<O ; l++){
                                if (B_matrix_local[i][l]==2){
                                    carbonyl[ca]=i;
                                    ca=ca+1;
                                    break;
                                }
                            }
                        }
                    }

                    // Construct residues from amine[] and carbonyl[]
                    for (int ai=0; ai<am; ai++){
                        int alphaC = amine[ai];
                        int aminoN = -1;
                        for (int n=C; n<N; n++){
                            int r = (alphaC < n) ? alphaC : n;
                            int c = (alphaC < n) ? n : alphaC;
                            if (B_matrix_local[r][c] == 1 && B_matrix_local[n][H] == 1){
                                aminoN = n;
                                break;
                            }                        
                        }
                        if (aminoN == -1) continue;

                        for (int ci=0; ci<ca; ci++){
                            int carC = carbonyl[ci];
                            int r2 = (alphaC < carC) ? alphaC : carC;
                            int c2 = (alphaC < carC) ? carC : alphaC;
                            if (B_matrix_local[r2][c2] == 1 && B_matrix_local[alphaC][H] >= 1){
                                if (!usedPair[alphaC][carC]){   //if(usedPair[alphaC][carC] != 1){
                                    residues[residueCount].carboxyC = carC;
                                    residues[residueCount].alphaC = alphaC;
                                    residues[residueCount].aminoN = aminoN;
                                    residueCount++;
                                    usedPair[alphaC][carC] = 1;
                                    break;  //一つのアミノ基は一つのカルボキシ基としかペアにならないようにするため、ペアが見つかったらループを抜ける
                                }
                            }
                        }
                    }
                    // N末端リストの作成
                    int NterminusCount = 0;
                    struct Residue Nterminus_list[50];

                    for (int ai=0; ai<am; ai++){
                        int alphaC = amine[ai];
                        int Nterminus = -1;
                        for (int n=C; n<N; n++){
                            int r = (alphaC < n) ? alphaC : n;
                            int c = (alphaC < n) ? n : alphaC;                        
                            if (B_matrix_local[r][c] == 1 && B_matrix_local[n][H] == 2){
                                Nterminus = n;
                                break;
                            }
                        }
                        if (Nterminus == -1) continue;

                        for (int ci=0; ci<ca; ci++){
                            int carC = carbonyl[ci];
                            int r2 = (alphaC < carC) ? alphaC : carC;
                            int c2 = (alphaC < carC) ? carC : alphaC;
                            if (B_matrix_local[r2][c2] == 1 && B_matrix_local[alphaC][H] >= 1){
                                if (!usedPair[alphaC][carC]){   //if(usedPair[alphaC][carC] != 1){
                                    Nterminus_list[NterminusCount].carboxyC = carC;
                                    Nterminus_list[NterminusCount].alphaC = alphaC;
                                    Nterminus_list[NterminusCount].aminoN = Nterminus;
                                    NterminusCount++;
                                    usedPair[alphaC][carC] = 1;
                                }
                            }
                        }
                    }

                    //ペプチド判定
                    int peptideBond[100][100];
                    int DP[50];  //重合度
                    for (int i=0; i<100; i++) for (int j=0; j<100; j++) peptideBond[i][j]=0;
                    for(int i=0; i<50; i++) DP[i]=0;

                    for(int i=0; i<NterminusCount; i++){
                        int Nterm_carC = Nterminus_list[i].carboxyC;
                        for(int j=0; j<residueCount; j++){
                            int res_aminoN = residues[j].aminoN;
                            if (B_matrix_local[Nterm_carC][res_aminoN] == 1){  //N末端残基(i)と他の残基(j)のC-N結合をチェック
                                peptideBond[i][j] = 1;
                                DP[i]=DP[i]+1;  //i番目のN末端残基に結合している残基数をカウント

                                //N末端残基(i)と他の残基(j)の間のC-N結合を判定後、チェーンを追跡
                                int currentResidueIdx = j;
                                int visitedResidues[100];
                                int visitedCount = 0;
                                visitedResidues[visitedCount] = j;
                                visitedCount++;

                                int chainComplete = 0;  //COOHが見つかったフラグ

                                while (!chainComplete && visitedCount < residueCount) {
                                    int curr_carC = residues[currentResidueIdx].carboxyC;
                                    int foundNext = 0;

                                    //現在の残基がCOOH基を持つか確認
                                    for(int l=0; l<CA; l++){
                                        if(curr_carC == carboxy[l]){
                                            //C末端（COOH基）を持つ残基を見つけた
                                            chainComplete = 1;
                                            DP[i]=DP[i]+1;  //C末端をカウント
                                            break;
                                        }
                                    }

                                    //COOHが見つからない場合、別の残基(k)とのC-N結合を探す
                                    if (!chainComplete) {
                                        for(int k=0; k<residueCount; k++){
                                            int isVisited = 0;
                                            for(int v=0; v<visitedCount; v++){
                                                if(visitedResidues[v] == k) {
                                                    isVisited = 1;
                                                    break;
                                                }
                                            }

                                            if (!isVisited) {
                                                int res_carC = residues[currentResidueIdx].carboxyC;
                                                int res_aminoN = residues[k].aminoN;
                                                //現在の残基のC末端と次の残基(k)のN末端をチェック
                                                if (B_matrix_local[res_carC][res_aminoN] == 1){
                                                    peptideBond[currentResidueIdx][k] = 1;
                                                    DP[i]=DP[i]+1;  //チェーン上の残基をカウント

                                                    //次の残基を現在の残基として設定
                                                    currentResidueIdx = k;
                                                    visitedResidues[visitedCount] = k;
                                                    visitedCount++;
                                                    foundNext = 1;
                                                    break;
                                                }
                                            }
                                        }

                                        if (!foundNext) {
                                            //次の残基が見つからない場合、重合度を０にしてループを終了
                                            DP[i] = 0;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (int i=0; i<NterminusCount; i++){
                        if (DP[i] == 2){
                            if(J==PostUV-1){
                                peptide_dimer_abundance_UV=peptide_dimer_abundance_UV+1;
                            }else if((STEP!=PostUV) && (J==STEP-1)){
                                peptide_dimer_abundance_PostUV=peptide_dimer_abundance_PostUV+1;
                            }else{
                                peptide_dimer_abundance[J/interval]=peptide_dimer_abundance[J/interval]+1;
                            }
                            run=I;
                            step=J;
                        }else if (DP[i] == 3){
                            if(J==PostUV-1){
                                peptide_trimer_abundance_UV=peptide_trimer_abundance_UV+1;
                            }else if((STEP!=PostUV) && (J==STEP-1)){
                                peptide_trimer_abundance_PostUV=peptide_trimer_abundance_PostUV+1;
                            }else{
                                peptide_trimer_abundance[J/interval]=peptide_trimer_abundance[J/interval]+1;
                            }
                            run=I;
                            step=J;
                        }else if (DP[i] >= 4){
                            if(J==PostUV-1){
                                peptide_oligomer_abundance_UV=peptide_oligomer_abundance_UV+1;
                            }else if((STEP!=PostUV) && (J==STEP-1)){
                                peptide_oligomer_abundance_PostUV=peptide_oligomer_abundance_PostUV+1;
                            }else{
                                peptide_oligomer_abundance[J/interval]=peptide_oligomer_abundance[J/interval]+1;
                            }
                            run=I;
                            step=J;
                        }
                    }
                }

                /*/ペプチド結合(アミド結合)  -CO-NH-
                if(J%20==0 || J>=PostUV-1){
                    for(int i=0;i<C;i=i+1){
                        for (int j=C ; j<N ; j++){
                            if (B_matrix_local[i][j]==1 && B_matrix_local[j][H]==1){//C-NH
                                for (int l=N ; l<O ; l++){
                                    if (B_matrix_local[i][l]==2){//C=O
                                        peptidelike_abundance[J]=peptidelike_abundance[J]+1;
                                        //run=I;
                                        //step=J;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }*/


                //Sugar
                /*if(J%20==0 || J>=PostUV){

                    //まずは行列内に存在する全ての分子の化学式を調べる
                    int aa=0;
                    for(int i=0;i<O;i=i+1){
                        for(int j=i+1; j<O; j=j+1){
                            if(B_matrix_local[i][j]!=0){//ノンゼロリスト
                                B_nz_local[aa].row1=i;
                                B_nz_local[aa].col1=j;
                                aa=aa+1;
                            }
                        }
                    }
                    int G[H+1];  //チェック未→０　済→１
                    int g[C][O];  //グループ（カッコ１）に原子のナンバー（カッコ２）を記録
                    int group=0,Total_atom[O];
                    int H_num[20],C_num[20],O_num[20],N_num[20];
                    for (int i=0 ; i<=H ; i++){
                        G[i]=0;
                    }

                    for (int i=0 ; i<aa ; i++){

                        if (B_nz_local[i].row1<C){

                            if(G[B_nz_local[i].row1]==0){

                                int atom = 0;
                                int old_atom = atom;

                                g[group][atom]=B_nz_local[i].row1; //グループに追加
                                G[B_nz_local[i].row1]=1; //追加済チェック
                                atom=atom+1;

                                g[group][atom]=B_nz_local[i].col1;
                                G[B_nz_local[i].col1]=1;
                                atom=atom+1;

                                B_nz_local[i].row1=H; // 以降無効化
                                B_nz_local[i].col1=H;

                                while (old_atom < atom){

                                    old_atom = atom ;

                                    for (int j=0 ; j<aa ; j++){

                                        if (G[B_nz_local[j].row1]==1 && G[B_nz_local[j].col1]==0){

                                            g[group][atom]=B_nz_local[j].col1;
                                            G[B_nz_local[j].col1]=1;
                                            atom=atom+1;

                                            B_nz_local[j].row1=H;
                                            B_nz_local[j].col1=H;

                                    }else if (G[B_nz_local[j].row1]==0 && G[B_nz_local[j].col1]==1){
                                        
                                            g[group][atom]=B_nz_local[j].row1;
                                            G[B_nz_local[j].row1]=1;
                                            atom=atom+1;

                                            B_nz_local[j].row1=H;
                                            B_nz_local[j].col1=H;

                                        }else if (G[B_nz_local[j].row1]==1 && G[B_nz_local[j].col1]==1){

                                            B_nz_local[j].row1=H;
                                            B_nz_local[j].col1=H;

                                        }
                                    }
                                }

                                Total_atom[group]=atom;
                                group = group + 1;
                            }
                        }
                    }

                    for (int i=0 ; i < group ; i++){

                        H_num[i]=0;
                        C_num[i]=0;
                        N_num[i]=0;
                        O_num[i]=0;

                        for (int j=0  ; j < Total_atom[i] ; j++){

                            H_num[i] = H_num[i] + B_matrix_local[g[i][j]][H] ;

                            if (0<=g[i][j] && g[i][j]<C){
                                C_num[i] = C_num[i] + 1 ;
                            }else if (C<=g[i][j] && g[i][j]<N){
                                N_num[i] = N_num[i] + 1 ;
                            }else if (N<=g[i][j] && g[i][j]<O){
                                O_num[i] = O_num[i] + 1 ;
                            }
                        }
                        //ここまででgroup番号iの分子の化学式が得られた

                        //以下、糖の判定
                        int CCbond=0,xx=0,oo=0;
                        if(C_num[i]>=3){  //C3以上のアルドース
                            if(C_num[i]==O_num[i] && 2*C_num[i]==H_num[i]){ //アルドースの判定(C×2n H×n O×2n)
                                if(N_num[i]==0){
                                    for (int k=0 ; k<Total_atom[i] ; k++){
                                        if(0<=g[i][k] && g[i][k]<C){ //i=group, k=atom
                                            for (int l=0 ; l<C ; l++){
                                                if(g[i][k]<l && B_matrix_local[g[i][k]][l]!=0){
                                                    CCbond = CCbond + 1;
                                                }else if(l<g[i][k] && B_matrix_local[l][g[i][k]]!=0){
                                                    CCbond = CCbond + 1;
                                                }
                                            }
                                        }
                                    }
                                    if(CCbond==2*(C_num[i]-1)){ //全てのC同士がつながっている
                                        for (int k=0 ; k<Total_atom[i] ; k++){
                                            if(0<=g[i][k] && g[i][k]<C){ //分子中の炭素kについて (g[i][k]=分子iを構成する原子k)
                                                int C_C=0, C_O=0, O_H=0;
                                                for (int l=0 ; l<C ; l++){ //各炭素原子が持つ他の炭素原子との結合次数を合計
                                                    if(g[i][k]<l){
                                                        C_C = C_C + B_matrix_local[g[i][k]][l];
                                                    }else if(l<g[i][k]){
                                                        C_C = C_C + B_matrix_local[l][g[i][k]];
                                                    }
                                            }
                                                for (int l=N ; l<O ; l++){ //各炭素原子が持つ酸素原子との結合次数を合計
                                                    C_O = C_O + B_matrix_local[g[i][k]][l];
                                                    if(B_matrix_local[g[i][k]][l]!=0){ //炭素につながる-OH基の数を合計
                                                        O_H = O_H + B_matrix_local[l][H];
                                                    }
                                                }
                                                if(C_C>2){ //C結合枝分かれ排除 (=炭素達が一本鎖でつながっている)
                                                    xx=1;
                                                    break;
                                                }else if(C_C==2){ //C-C-C <-真ん中の炭素について
                                                    if(C_O!=1 || O_H!=1){ //Oがつながってない炭素、C=O、C-O
                                                        xx=1;
                                                        break;
                                                    }
                                                }else if(C_C==1){ //(環状)C-CH2(OH) または C-CH(OH)-O, (鎖状)C-CH=O または C-CH2(OH)
                                                    if(C_O==0){ //(環状)C-CH2(OH) または C-CH(OH)-O, (鎖状)C-CH=O または C-CH2(OH)
                                                        xx=1;
                                                        break;
                                                    }else if(C_O==2){ //C-O-O-Cを防ぐため (少なくとも一つはC=O結合持っている)
                                                        oo=1;
                                                    }
                                                }
                                            }
                                        }
                                        if(xx==0 && oo==1){
                                            sugar_abundance[J]=sugar_abundance[J]+1; // 鎖状
                                            //aa=0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }*/
                //--------------------------------分子種の判別おわり--------------------------------



                //--------------------------------分子描画(Marvin view)--------------------------------
                if(I==0){
                    int aa=0;
                    int Total_1=0, Total_2=0;
                    //ラジカル以外
                    for(int i=0;i<S;i=i+1){
                        for(int j=i+1; j<S; j=j+1){
                            if(B_matrix_local[i][j]!=0){
                                B_nz_local[aa].row1=i;
                                B_nz_local[aa].col1=j;
                                aa=aa+1;
                            }
                        }
                    }
                    Total_1=aa;
                    //ラジカル(R-X)
                    for(int i=0;i<S;i=i+1){
                        for(int j=0; j<B_matrix_local[i][X]; j=j+1){
                            B_nz_local[aa].row1=i;
                            aa=aa+1;
                        }
                    }
                    Total_2=aa;


                    //以下、Marvin View用の出力ファイル
                    fprintf(s,"\n");

                    fprintf(s,"  Yoko   Run:%d Step:%d\n",I+1, J+1);

                    fprintf(s,"\n");

                    fprintf(s,"%3d%3d  0  0  0  0            999 V2000\n",S+Total_2-Total_1+B_matrix_local[H][X]+B_matrix_local[X][X]*2,Total_2+B_matrix_local[X][X]); //H以外の原子の数，H-R以外の結合の総数

                    for (int i=0 ; i<C ;i++){
                        fprintf(s,"    0.0       0.0       0.0    C   0  0  0  0  0  0  0  0  0  0  0  0\n");
                    }
                    for (int i=C ; i<N ; i++){
                        fprintf(s,"    0.0       0.0       0.0    N   0  0  0  0  0  0  0  0  0  0  0  0\n");
                    }
                    for (int i=N ; i<O ; i++){
                        fprintf(s,"    0.0       0.0       0.0    O   0  0  0  0  0  0  0  0  0  0  0  0\n");
                    }
                    for (int i=O ; i<S ; i++){
                        fprintf(s,"    0.0       0.0       0.0    S   0  0  0  0  0  0  0  0  0  0  0  0\n");
                    }
                    for (int i=0 ; i<Total_2-Total_1+B_matrix_local[H][X]+B_matrix_local[X][X]*2 ; i++){
                        fprintf(s,"    0.0       0.0       0.0    F   0  0  0  0  0  0  0  0  0  0  0  0\n");
                    }
                    for (int i=0 ; i<Total_1 ; i++){
                        fprintf(s,"%3d%3d%3d  0  0  0  0\n",B_nz_local[i].row1+1,B_nz_local[i].col1+1,B_matrix_local[B_nz_local[i].row1][B_nz_local[i].col1]);
                        /////////①原子の番号１個目、②２個目、③結合次数
                    }
                    for (int i=Total_1 ; i<Total_2 ; i++){
                        fprintf(s,"%3d%3d%3d  0  0  0  0\n",B_nz_local[i].row1+1,S+1+i-Total_1,1); //R-X
                } 
                    for (int i=0 ; i<B_matrix_local[X][X] ; i++){
                        fprintf(s,"%3d%3d%3d  0  0  0  0\n",S+Total_2-Total_1+1+2*i,S+Total_2-Total_1+1+2*i+1,1); //X-X
                    } 

                    fprintf(s,"M  END\n$$$$\n");

                    abc=abc+1;
                }
                //--------------------------------分子描画おわり--------------------------------*/

            //#pragma omp critical
            //printf("Run:%6d Step:%10d\r",I,J);
            //fflush(stdout);

            } //　J : ステップ

            // ★ メモリ解放
                free(B_matrix_local);
                free(B_nz_local);
                free(R_nz_local);
                free(D_energy_local);
                free(dH_local);
                free(dS_local);
                free(P_local);

        } //　I : ラン
    } // #pragma omp parallel 終了

    fprintf(f,"%d %lf %d %d %lf %lf %lf %lf\n", 0, 0.0, 0, 0, (Gly+Ala+Ser+Cys)*1.0, 0.0, 0.0, 0.0);

    for (int i=0 ; i<PostUV/interval ; i=i+1){
        if(i==PostUV-1) continue;
        fprintf(f,"%d %lf %d %d %lf %lf %lf %lf\n",i*interval+1, (i*interval+1)*dt, collision_ch[i], collision_th[i], aacid_abundance[i]/(double)TRIAL, peptide_dimer_abundance[i]/(double)TRIAL, peptide_trimer_abundance[i]/(double)TRIAL, peptide_oligomer_abundance[i]/(double)TRIAL);
    }
    fprintf(f,"%d %lf %d %d %lf %lf %lf %lf\n",PostUV, (PostUV)*dt, collision_ch_UV, collision_th_UV, aacid_abundance_UV/(double)TRIAL, peptide_dimer_abundance_UV/(double)TRIAL, peptide_trimer_abundance_UV/(double)TRIAL, peptide_oligomer_abundance_UV/(double)TRIAL);
    for (int i=PostUV/interval ; i<STEP/interval ; i=i+1){
        if(i==STEP-1) continue;
        fprintf(f,"%d %lf %d %d %lf %lf %lf %lf\n",i*interval+1, dt*PostUV+(i*interval+1-PostUV)*dt, collision_ch[i], collision_th[i], aacid_abundance[i]/(double)TRIAL, peptide_dimer_abundance[i]/(double)TRIAL, peptide_trimer_abundance[i]/(double)TRIAL, peptide_oligomer_abundance[i]/(double)TRIAL);
    }
    if(STEP!=PostUV) fprintf(f,"%d %lf %d %d %lf %lf %lf %lf\n",STEP, (STEP)*dt, collision_ch_PostUV, collision_th_PostUV, aacid_abundance_PostUV/(double)TRIAL, peptide_dimer_abundance_PostUV/(double)TRIAL, peptide_trimer_abundance_PostUV/(double)TRIAL, peptide_oligomer_abundance_PostUV/(double)TRIAL);

    /*/計算時間 出力ファイル
    FILE* t = fopen("file_time_10Gly_test_calculation.txt", "w");

    long cpu_time;
    double seconds;

    cpu_time = clock();
    seconds = (double)cpu_time / CLOCKS_PER_SEC;    //秒に直す

    printf("CPU Time: %ld ticks\n", cpu_time);
    printf("Elapsed Time: %f seconds\n", seconds);

    fprintf(t, "CPU Time: %ld ticks\n", cpu_time);
    fprintf(t, "Elapsed Time: %f seconds\n", seconds);
    fprintf(t, "              %f minutes\n", seconds / 60.0);
    fprintf(t, "              %f hours\n", seconds / 3600.0);

    fclose(t);//*/
    fclose(f);
    fclose(s);

    printf("All calculations have been completed!\n");

    return 0;
}
                