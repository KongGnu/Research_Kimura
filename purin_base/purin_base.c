#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

// SDF出力は非常に重いためデフォルト無効。必要なら1に変更。
#ifndef WRITE_SDF
#define WRITE_SDF 0
#endif
 
int C,N,O,H,X;

//Initial material （他の分子を使いたい場合は以下のコードを書き換えてください）
int CH3OH = 2;
int CH2O = 5;
int HCOOH = 0;
int HCN = 0;
int NH3 = 9;
int H2O = 22;
int X2 = 1;
int PURINE = 3;  // 初期プリン分子 C5H4N4 の個数

void COUNT(int a ,int b, int c, int d, int e, int f, int g){ //炭素原子の行（列）＝0行目〜C-1行目、窒素分子の行（列）＝C行目〜N-1行目、酸素原子の行（列）＝N行目〜O-1行目、水素原子の行（列）＝O行目=H行目、X原子の行（列）＝H+1行目
    C=0,O=0,H=0,N=0;

    //Carbon
    C = a + b + c + d + 5*PURINE;

    //Nitrogen
    N = C + d + e + 4*PURINE;

    //Oxygen
    O = N + a + b + c*2 + f;
    
    //Hydrogen
    H = O ;

    X = H + 1;
}

void INIT(int a[80][80]){ //初期分子の行列設定

    COUNT(CH3OH,CH2O,HCOOH,HCN,NH3,H2O,X2);

    O = N ;
    N = C ;
    C = 0 ;

    for(int i=0 ; i<=X ; i++){
        for(int j=i ; j<=X ; j++){
            a[i][j]=0;
        }
    }



    //CH3OH
    for (int i=0 ; i<CH3OH ; i++){
        a[C][H]=3;
        a[C][O]=1;
        a[O][H]=1;
        C=C+1;
        O=O+1;
    }

    //CH2O
    for (int i=0 ; i<CH2O ; i++){
        a[C][O]=2;
        a[C][H]=2;
        C=C+1;
        O=O+1;
    }

    //HCOOH
    for (int i=0 ; i<HCOOH ; i++){
        a[C][O]=2;
        a[C][H]=1;
        a[C][O+1]=1;
        a[O+1][H]=1;
        C=C+1;
        O=O+2;
    }

    //HCN
    for (int i=0 ; i<HCN ; i++){
        a[C][N]=3;
        a[C][H]=1;
        C=C+1;
        N=N+1;
    }

    //NH3
    for (int i=0 ; i<NH3 ; i++){
        a[N][H]=3;
        N=N+1;
    }

    // プリン C5H4N4：6員環+5員環の縮合骨格（Kekule 型）
    // 骨格番号: 0,1,2,4,7 は C / 3,5,6,8 は N
    // 環: 0-1-2-3-4-5-0 および 0-1-6-7-8-0
    for (int p=0; p<PURINE; p++) {
        int atom[9] = {C, C+1, C+2, N, C+3, N+1, N+2, C+4, N+3};
        const int edges[10][3] = {
            {0,1,1}, {1,2,2}, {2,3,1}, {3,4,2}, {4,5,1},
            {5,0,2}, {1,6,1}, {6,7,2}, {7,8,1}, {8,0,1}
        };
        for (int k=0; k<10; k++) {
            int u=atom[edges[k][0]], v=atom[edges[k][1]];
            a[u<v ? u:v][u<v ? v:u] = edges[k][2];
        }
        // C2, C4, C8 にそれぞれ H; N9 に H
        a[atom[2]][H]=1;
        a[atom[4]][H]=1;
        a[atom[7]][H]=1;
        a[atom[8]][H]=1;
        C+=5;
        N+=4;
    }

    //H2O
    for (int i=0 ; i<H2O ; i++ ){
        a[O][H]=2;
        O=O+1;
    }

    //X2
    a[X][X]=X2;

    //H2
    a[H][H]=0;
}

void ENERGY(double b[50][80][80]){ //結合エネルギー

    COUNT(CH3OH,CH2O,HCOOH,HCN,NH3,H2O,X2);

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

    //O-H結合
    for (int i=N ; i<O ; i++){
        b[1][i][H]=464.;
        b[2][i][H]=464.*2.;
        b[3][i][H]=464.*3.;

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
        b[i][X][X]=-1000.*(double)i;
    }
}

// プリン骨格（縮合6員環+5員環）の9原子だけを検出する。
// 6員環: 0-1-2-3-4-5-0  (C-C-C-N-C-N)
// 5員環: 0-1-6-7-8-0  (C-C-N-C-N)
// 0,1は共有される2つの炭素。環内の原子総数 C5N4。
// H, O, 置換基、および分子全体の化学式には制限を課さない。
// 同じ分子内で複数の対応付けが見つかっても1回だけカウントする。
static const int purine_types[9] = {0,0,0,1,0,1,1,0,1}; // 0=C, 1=N
static const int purine_edges[10][2] = {
    {0,1},{1,2},{2,3},{3,4},{4,5},{5,0},
    {1,6},{6,7},{7,8},{8,0}
};

static int heavy_bond(const int b[80][80], int x, int y) {
    return x < y ? b[x][y] : b[y][x];
}

static int purine_edge_expected(int u, int v) {
    for (int e=0; e<10; e++) {
        if ((purine_edges[e][0]==u && purine_edges[e][1]==v) ||
            (purine_edges[e][0]==v && purine_edges[e][1]==u)) return 1;
    }
    return 0;
}

static int match_purine_core(const int b[80][80], const int atoms[80],
                             int count, int depth, int selected[9],
                             unsigned char used[80]) {
    if (depth == 9) return 1;
    for (int i=0; i<count; i++) {
        int a=atoms[i];
        if (a>=N || used[a]) continue; // C/N のみ環原子として許可
        if ((a>=C) != purine_types[depth]) continue;
        int valid=1;
        for (int j=0; j<depth; j++) {
            // 環骨格の内部に余分な結合（chord）があるものを除外
            if ((heavy_bond(b,a,selected[j])>0) !=
                purine_edge_expected(depth,j)) {
                valid=0; break;
            }
        }
        if (!valid) continue;
        used[a]=1;
        selected[depth]=a;
        if (match_purine_core(b,atoms,count,depth+1,selected,used)) return 1;
        used[a]=0;
    }
    return 0;
}

static int detect_purine_cores(const int b[80][80]) {
    unsigned char visited[80]={0};
    int count_cores=0;
    // 重原子結合で連結した「分子」単位で判定する
    for (int start=0; start<O; start++) {
        if (visited[start]) continue;
        int atoms[80], n=0, head=0;
        atoms[n++]=start;
        visited[start]=1;
        while (head<n) {
            int v=atoms[head++];
            for (int w=0; w<O; w++) {
                if (!visited[w] && w!=v && heavy_bond(b,v,w)>0) {
                    visited[w]=1;
                    atoms[n++]=w;
                }
            }
        }
        int nc=0, nn=0;
        for (int i=0; i<n; i++) {
            nc += atoms[i]<C;
            nn += atoms[i]>=C && atoms[i]<N;
        }
        if (nc<5 || nn<4) continue;
        int selected[9]={0};
        unsigned char used[80]={0};
        if (match_purine_core(b,atoms,n,0,selected,used)) count_cores++;
    }
    return count_cores;
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


int main() {
    printf("start\n");
    fflush(stdout);
    FILE* f = fopen("file_abundance.txt", "w");//アバンダンス出力ファイル
    FILE* s = WRITE_SDF ? fopen("file_marvinview.sdf", "w") : NULL;
    if (!f || (WRITE_SDF && !s)) { perror("fopen"); return 1; }

    //Parameters
    int TRIAL=10000;    //試行回数
    double T=100.;      //温度 K
    int STEP=1000;      //全ステップ数
    int PostUV=900;    //PostUVphase開始ステップ
    double alpha=1.;
    double beta=100.; //[kJ/mol]
    double timescale=1000.; //yr
    double Ea=40*(T/100)*(1.1+0.048*log10((timescale/1000.)*(T/100)));//[kJ/mol]：ある温度Tにおいて、timescaleの期間内に起きる反応の活性化エネルギーの最大値

    // INITで原子の境界 C,N,O,H,X を確定してから並列領域へ
    int initial_matrix[80][80] = {{0}};
    INIT(initial_matrix);
    double (*B_energy)[80][80] = calloc(50, sizeof(*B_energy));
    if (!B_energy) { perror("calloc"); fclose(f); if(s) fclose(s); return 1; }
    if (X>=80) {
        fprintf(stderr, "ERROR: too many atoms for B_matrix[80][80] (X=%d)\n", X);
        fclose(f); if(s) fclose(s);
        return 1;
    }
    ENERGY(B_energy);

    // スレッドに依存しないtrial別乱数シード
    const double ra=48271.;
    const double rm=2147483647.;

    //気体定数
    double R=0.00831; 

    //Abundance
    double aacid_abundance[5000];
    double sugar_abundance[5000];
    double purine_core_abundance[5000] = {0};
    for (int i=0 ; i<STEP ; i++){
        aacid_abundance[i]=0;
        sugar_abundance[i]=0;
    }

    int initial_cores = detect_purine_cores(initial_matrix);
    printf("initial purine cores = %d (expected %d)\n", initial_cores, PURINE);
    if (initial_cores != PURINE) {
        fprintf(stderr, "ERROR: initial purine core detection mismatch\n");
        fclose(f); if(s) fclose(s); free(B_energy); return 1;
    }

    // 各スレッドが作業配列を専有（大きい配列をスレッドスタックに置かない）
    #pragma omp parallel
    {
        int B_matrix[80][80];
        struct nonzero *B_nz = malloc(10000 * sizeof(*B_nz));
        struct nonzero *R_nz = malloc(100000 * sizeof(*R_nz));
        struct nonzero CO_dou[100], CO_tri[100];
        double *D_energy = malloc(100000 * sizeof(*D_energy));
        double *dH = malloc(100000 * sizeof(*dH));
        if (!B_nz || !R_nz || !D_energy || !dH) {
            fprintf(stderr, "ERROR: worker allocation failed\n");
            abort();
        }
        #pragma omp for schedule(dynamic)
        for (int I=0 ; I<TRIAL ; I++) {
            memcpy(B_matrix, initial_matrix, sizeof(B_matrix));
            double r0 = (double)(1ULL + ((7788776655ULL + 104729ULL*(uint64_t)I) % 2147483646ULL));

        for (int J=0 ; J<STEP; J++){ //ステップ回数

            if(J<PostUV){ 
                B_matrix[X][X]=X2; //UVphaseの間はX2を含める
            }else{
                B_matrix[X][X]=0; //PostUVphaseはX2なし
            }

            // B_nzはスレッド専用のヒープ配列

            // CO_dou, CO_triもスレッド専用

            int Total_nz=0; //ノンゼロ成分の個数
            int CO2=0, CO3=0; //CO二重結合, CO三重結合の個数

            //ノンゼロ成分の行(row1)、列(col1)を記録
            for(int i=0;i<=H;i=i+1){ //C~H列の探索(X2はリストに含めない)
                for(int j=i; j<=X; j=j+1){ //C~X行の探索
                    if(j==H || j==X){ 
                        for (int k=0 ; k<B_matrix[i][j] ; k++){ 
                            B_nz[Total_nz].row1=i;
                            B_nz[Total_nz].col1=j;
                            Total_nz=Total_nz+1;
                        }
                    }else if (i<C && j>=N){ 
                        if(B_matrix[i][j]==1 || B_matrix[i][j]==2){//CO分子はカウントしないように
                            B_nz[Total_nz].row1=i;
                            B_nz[Total_nz].col1=j;
                            Total_nz=Total_nz+1;
                        }
                    }else if(B_matrix[i][j]!=0){
                        B_nz[Total_nz].row1=i;
                        B_nz[Total_nz].col1=j;
                        Total_nz=Total_nz+1;
                    }
                }
            }

            //CO二重結合、三重結合成分の行(row1)、列(col1)を記録
            for(int i=0;i<C;i=i+1){
                for(int j=N; j<O; j=j+1){
                    if (B_matrix[i][j]==2){
                        CO_dou[CO2].row1=i;
                        CO_dou[CO2].col1=j;
                        CO2=CO2+1;
                    }else if (B_matrix[i][j]==3){
                        CO_tri[CO3].row1=i;
                        CO_tri[CO3].col1=j;
                        CO3=CO3+1;
                    }
                }
            }


            // R_nzはスレッド専用のヒープ配列
            int bb=0; //R行列の数＝反応候補の数
            //R_nz[bb].row1, R_nz[bb].col1: 切断する結合の行、列番号
            //R_nz[bb].row2, R_nz[bb].col2: 切断するもう一つの結合の行、列番号
            //R_nz[bb].row3, R_nz[bb].col3: 生成する結合の行、列番号
            //R_nz[bb].row4, R_nz[bb].col4: 生成するもう一つの結合の行、列番号

            // D_energy, dHはスレッド専用のヒープ配列


            //変化無 NO CHANGE
            R_nz[bb].row1=H, R_nz[bb].col1=H; 
            R_nz[bb].row2=H, R_nz[bb].col2=H;
            R_nz[bb].row3=H, R_nz[bb].col3=H;
            R_nz[bb].row4=H, R_nz[bb].col4=H;
            D_energy[bb]=Ea;
            bb=bb+1;

            //結合行列のノンゼロリストをもとに、切断させる結合・生成させる結合の組み合わせを作成（反応候補リスト作成）
            for (int l=0 ; l<Total_nz ; l++){
                for (int m=l+1 ; m<Total_nz ; m++){ //l番目以降のm番目のノンゼロ成分について
                    if(B_nz[l].row1!=B_nz[m].row1){  //同じ原子が持つ結合２本は選べない
                        if(B_nz[l].col1==H || B_nz[l].col1==X){  //R-HまたはR-Xの場合B_nz[l].col1とB_nz[m].col1は同じでもいい

                            R_nz[bb].row1=B_nz[l].row1, R_nz[bb].col1=B_nz[l].col1;
                            R_nz[bb].row2=B_nz[m].row1, R_nz[bb].col2=B_nz[m].col1;
                            R_nz[bb].row3=B_nz[l].row1, R_nz[bb].col3=B_nz[m].row1;
                            if(B_nz[l].col1<=B_nz[m].col1){
                                R_nz[bb].row4=B_nz[l].col1, R_nz[bb].col4=B_nz[m].col1;
                            }else{
                                R_nz[bb].row4=B_nz[m].col1, R_nz[bb].col4=B_nz[l].col1;
                            }
                            bb = bb + 1;

                            R_nz[bb].row1=B_nz[l].row1, R_nz[bb].col1=B_nz[l].col1;
                            R_nz[bb].row2=B_nz[m].row1, R_nz[bb].col2=B_nz[m].col1;
                            R_nz[bb].row3=B_nz[l].row1, R_nz[bb].col3=B_nz[m].col1;
                            if(B_nz[m].row1<=B_nz[l].col1){
                                R_nz[bb].row4=B_nz[m].row1, R_nz[bb].col4=B_nz[l].col1;
                            }else{
                                R_nz[bb].row4=B_nz[l].col1, R_nz[bb].col4=B_nz[m].row1;
                            }
                            bb = bb + 1;

                        }else if(B_nz[l].col1!=B_nz[m].row1 && B_nz[l].col1!=B_nz[m].col1){ //R1-R2 + R3-R4ならR2はR3,R4と異なる原子でなくてはならない
                            if(B_nz[l].col1<B_nz[m].col1){
                        
                                R_nz[bb].row1=B_nz[l].row1, R_nz[bb].col1=B_nz[l].col1;
                                R_nz[bb].row2=B_nz[m].row1, R_nz[bb].col2=B_nz[m].col1;
                                R_nz[bb].row3=B_nz[l].row1, R_nz[bb].col3=B_nz[m].row1;
                                R_nz[bb].row4=B_nz[l].col1, R_nz[bb].col4=B_nz[m].col1;

                                bb = bb + 1;

                            }else{

                                R_nz[bb].row1=B_nz[l].row1, R_nz[bb].col1=B_nz[l].col1;
                                R_nz[bb].row2=B_nz[m].row1, R_nz[bb].col2=B_nz[m].col1;
                                R_nz[bb].row3=B_nz[l].row1, R_nz[bb].col3=B_nz[m].row1;
                                R_nz[bb].row4=B_nz[m].col1, R_nz[bb].col4=B_nz[l].col1;

                                bb = bb + 1;
                                
                            }

                            if(B_nz[l].col1<B_nz[m].row1){
                    
                                R_nz[bb].row1=B_nz[l].row1, R_nz[bb].col1=B_nz[l].col1;
                                R_nz[bb].row2=B_nz[m].row1, R_nz[bb].col2=B_nz[m].col1;
                                R_nz[bb].row3=B_nz[l].row1, R_nz[bb].col3=B_nz[m].col1;
                                R_nz[bb].row4=B_nz[l].col1, R_nz[bb].col4=B_nz[m].row1;

                                bb = bb + 1;

                            }else{

                                R_nz[bb].row1=B_nz[l].row1, R_nz[bb].col1=B_nz[l].col1;
                                R_nz[bb].row2=B_nz[m].row1, R_nz[bb].col2=B_nz[m].col1;
                                R_nz[bb].row3=B_nz[l].row1, R_nz[bb].col3=B_nz[m].col1;
                                R_nz[bb].row4=B_nz[m].row1, R_nz[bb].col4=B_nz[l].col1;

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
                        if (B_matrix[k][CO_dou[i].row1]==2){
                            break;
                        }else if (B_matrix[k][CO_dou[i].row1]==1){
                            for (int l=k+1 ; l<H ; l++){
                                if (CO_dou[i].row1<l && B_matrix[CO_dou[i].row1][l]==1){
                                    R_nz[bb].row1=k             , R_nz[bb].col1=CO_dou[i].row1;
                                    R_nz[bb].row2=CO_dou[i].row1, R_nz[bb].col2=l;
                                    R_nz[bb].row3=k             , R_nz[bb].col3=l;
                                    R_nz[bb].row4=CO_dou[i].row1, R_nz[bb].col4=CO_dou[i].col1;
                                    bb=bb+1;
                                    break;
                                }else if (l<CO_dou[i].row1 && B_matrix[l][CO_dou[i].row1]==1){
                                    R_nz[bb].row1=k             , R_nz[bb].col1=CO_dou[i].row1;
                                    R_nz[bb].row2=l             , R_nz[bb].col2=CO_dou[i].row1;
                                    R_nz[bb].row3=k             , R_nz[bb].col3=l;
                                    R_nz[bb].row4=CO_dou[i].row1, R_nz[bb].col4=CO_dou[i].col1;
                                    bb=bb+1;
                                    break;
                                }
                            }
                        }
                    }else if (CO_dou[i].row1<k){
                        if (B_matrix[CO_dou[i].row1][k]==2){
                            break;
                        }else if (B_matrix[CO_dou[i].row1][k]==1){
                            for (int l=k+1 ; l<H ; l++){
                                if (B_matrix[CO_dou[i].row1][l]==1){
                                    R_nz[bb].row1=CO_dou[i].row1, R_nz[bb].col1=k;
                                    R_nz[bb].row2=CO_dou[i].row1, R_nz[bb].col2=l;
                                    R_nz[bb].row3=k             , R_nz[bb].col3=l;
                                    R_nz[bb].row4=CO_dou[i].row1, R_nz[bb].col4=CO_dou[i].col1;
                                    bb=bb+1;
                                    break;
                                }
                            }
                        }
                    }
                    for(int k=H ; k<=X ; k++){
                        if (B_matrix[CO_dou[i].row1][k]==2){
                            R_nz[bb].row1=CO_dou[i].row1, R_nz[bb].col1=k;
                            R_nz[bb].row2=CO_dou[i].row1, R_nz[bb].col2=k;
                            R_nz[bb].row3=k,              R_nz[bb].col3=k;
                            R_nz[bb].row4=CO_dou[i].row1, R_nz[bb].col4=CO_dou[i].col1;
                            bb=bb+1;
                            break;
                        }
                    }
                }
            }

            //CO分解反応 CO + R-R' -> R-CO-R'
            for(int i=0 ; i<CO3 ; i++){ 
                for (int l=0 ; l<Total_nz ; l++){
                    R_nz[bb].row1=B_nz[l].row1,   R_nz[bb].col1=B_nz[l].col1;
                    R_nz[bb].row2=CO_tri[i].row1, R_nz[bb].col2=CO_tri[i].col1;

                    if(CO_tri[i].row1<B_nz[l].row1){

                        R_nz[bb].row3=CO_tri[i].row1, R_nz[bb].col3=B_nz[l].row1;
                        R_nz[bb].row4=CO_tri[i].row1, R_nz[bb].col4=B_nz[l].col1;
                        bb = bb + 1;

                    }else if (B_nz[l].row1<CO_tri[i].row1 && CO_tri[i].row1<B_nz[l].col1){

                        R_nz[bb].row3=B_nz[l].row1,   R_nz[bb].col3=CO_tri[i].row1;
                        R_nz[bb].row4=CO_tri[i].row1, R_nz[bb].col4=B_nz[l].col1;
                        bb = bb + 1;

                    }else if (B_nz[l].col1<CO_tri[i].row1){

                        R_nz[bb].row3=B_nz[l].row1,   R_nz[bb].col3=CO_tri[i].row1;
                        R_nz[bb].row4=B_nz[l].col1,   R_nz[bb].col4=CO_tri[i].row1;
                        bb = bb + 1;
                    }
                }
                for(int i=0 ; i<CO3 ; i++){ //ラジカルとの反応
                    for (int j=0 ; j<X ; j++){
                        for (int k=0 ; k<B_matrix[j][X] ; k++){
                            R_nz[bb].row1=j             , R_nz[bb].col1=X;
                            R_nz[bb].row2=CO_tri[i].row1, R_nz[bb].col2=CO_tri[i].col1;
                            R_nz[bb].row4=CO_tri[i].row1, R_nz[bb].col4=X;
                            if(CO_tri[i].row1<j){
                                R_nz[bb].row3=CO_tri[i].row1, R_nz[bb].col3=j;
                                bb = bb + 1;
                            }else{
                                R_nz[bb].row3=j, R_nz[bb].col3=CO_tri[i].row1;
                                bb = bb + 1;
                            }
                        }
                    }
                }
            }

            //光解離反応 R-R' + X2 -> R-X + R'-X
            for (int l=0 ; l<Total_nz ; l++){
                for (int j=0 ; j<B_matrix[X][X] ; j++){
                    if(B_nz[l].col1!=X){
                        R_nz[bb].row1=B_nz[l].row1, R_nz[bb].col1=B_nz[l].col1;
                        R_nz[bb].row2=X,            R_nz[bb].col2=X;
                        R_nz[bb].row3=B_nz[l].row1, R_nz[bb].col3=X;
                        R_nz[bb].row4=B_nz[l].col1, R_nz[bb].col4=X;
                        D_energy[bb]=0;//光解離反応は活性化エネルギーなし
                        bb=bb+1;
                    }
                }
            }

            for (int i=1 ; i<bb ; i++){
                dH[i]=B_energy[B_matrix[R_nz[i].row1][R_nz[i].col1]][R_nz[i].row1][R_nz[i].col1]//切る
                            -B_energy[B_matrix[R_nz[i].row1][R_nz[i].col1]-1][R_nz[i].row1][R_nz[i].col1]

                            +B_energy[B_matrix[R_nz[i].row2][R_nz[i].col2]][R_nz[i].row2][R_nz[i].col2]//切る
                            -B_energy[B_matrix[R_nz[i].row2][R_nz[i].col2]-1][R_nz[i].row2][R_nz[i].col2]

                            +B_energy[B_matrix[R_nz[i].row3][R_nz[i].col3]][R_nz[i].row3][R_nz[i].col3]
                            -B_energy[B_matrix[R_nz[i].row3][R_nz[i].col3]+1][R_nz[i].row3][R_nz[i].col3]//結合

                            +B_energy[B_matrix[R_nz[i].row4][R_nz[i].col4]][R_nz[i].row4][R_nz[i].col4]
                            -B_energy[B_matrix[R_nz[i].row4][R_nz[i].col4]+1][R_nz[i].row4][R_nz[i].col4];//結合 

                if(R_nz[i].row1==X && R_nz[i].col1==X){//光解離反応
                    dH[i]=100.;//100は適当な値です
                }else if(R_nz[i].row2==X && R_nz[i].col2==X){//光解離反応
                    dH[i]=100.;//100は適当な値です
                }else if(R_nz[i].row4==X && R_nz[i].col4==X){//ラジカル-ラジカル反応
                    dH[i]=dH[i]+B_energy[1][X][X];
                    D_energy[i]=dH[i];
                }else if(R_nz[i].row3==X && R_nz[i].col3==X){//ラジカル-ラジカル反応
                    dH[i]=dH[i]+B_energy[1][X][X];
                    D_energy[i]=dH[i];
                }else{
                    D_energy[i]=alpha*dH[i]+beta;
                }
                
                if(D_energy[i]<0){//活性化エネルギーはマイナスにならないように
                    D_energy[i]=0;
                }
            }

            //HX+HX->H2+X2
            for(int i=0 ; i<(B_matrix[H][X]*(B_matrix[H][X]-1))/2 ; i++){
                R_nz[bb].row1=H, R_nz[bb].col1=X;
                R_nz[bb].row2=H, R_nz[bb].col2=X;
                R_nz[bb].row3=H, R_nz[bb].col3=H;
                R_nz[bb].row4=X, R_nz[bb].col4=X;
                D_energy[bb]=0;//activation E
                dH[bb]=-B_energy[1][H][H];
                bb=bb+1;
            }

            int Total_R1=bb;

            //C三O+XOH->HO-CX=O
            for(int i=0 ; i<CO3 ; i++){ 
                for (int j=N ; j<O ; j++){
                    if(B_matrix[j][H]==1 && B_matrix[j][X]==1){
                        R_nz[bb].row1=j             , R_nz[bb].col1=X;
                        R_nz[bb].row2=CO_tri[i].row1, R_nz[bb].col2=CO_tri[i].col1;
                        R_nz[bb].row3=CO_tri[i].row1, R_nz[bb].col3=j;
                        R_nz[bb].row4=CO_tri[i].row1, R_nz[bb].col4=X;
                        D_energy[bb]=100;//activation E
                        dH[bb]=-31.;
                        bb = bb + 1;
                    }
                }
            }

            int Total_R=bb;

            //min=活性化エネルギーの最小値
            double min = D_energy[0];
            for (int i=1 ; i<Total_R ; i++){
                if ( D_energy[i] < min ){
                    min = D_energy[i];
                }
                if(min==0.){
                    break;
                }
            }

            //反応速度の和
            double Total_exp=0;
            for (int i=0 ; i<Total_R ; i++){
                if ((D_energy[i]-min)<590.*(T/100.)){//underflow防止
                    Total_exp = Total_exp + exp(-(D_energy[i]-min)/(R*T));
                }
            }
            
            //ランダムに反応選択 → E行列（生成物行列）生成
            int re;
            double ex ;
            ex = 0 ;
            r0 = fmod(ra*r0,rm); //乱数
            for (int i=0 ; i<Total_R ; i++ ){
                if ((D_energy[i]-min)<590.*(T/100.)){
                    if (ex <= r0/rm && r0/rm < ex + exp(-(D_energy[i]-min)/(R*T))/Total_exp){
                        re=i;//選ばれた反応の番号
                        B_matrix[R_nz[i].row3][R_nz[i].col3] = B_matrix[R_nz[i].row3][R_nz[i].col3] + 1;//made
                        B_matrix[R_nz[i].row4][R_nz[i].col4] = B_matrix[R_nz[i].row4][R_nz[i].col4] + 1;//made
                        B_matrix[R_nz[i].row1][R_nz[i].col1] = B_matrix[R_nz[i].row1][R_nz[i].col1] - 1;//broken
                        B_matrix[R_nz[i].row2][R_nz[i].col2] = B_matrix[R_nz[i].row2][R_nz[i].col2] - 1;//broken
                        if(i>=Total_R1){
                            B_matrix[R_nz[i].row3][R_nz[i].col3] = B_matrix[R_nz[i].row3][R_nz[i].col3] + 1;//C=Omade
                            B_matrix[R_nz[i].row4][R_nz[i].col4] = B_matrix[R_nz[i].row4][R_nz[i].col4] - 1;//C-Xbroken
                            B_matrix[R_nz[i].row1][H]            = B_matrix[R_nz[i].row1][H]            - 1;//O-Hbroken
                            B_matrix[H][X]                       = B_matrix[H][X]                       + 1;//H-Xmade
                        }
                        break;
                    }
                    ex = ex + exp(-(D_energy[i]-min)/(R*T))/Total_exp;
                }
            }


            //--------------------------------分子種の判別--------------------------------
            if (J%1==0 || J>=PostUV) {
                int ncores = detect_purine_cores(B_matrix);
                #pragma omp atomic update
                purine_core_abundance[J] += ncores;
            }

            #if WRITE_SDF
            // SDFのレコード全体を排他出力してファイル破損を防ぐ
            #pragma omp critical(sdf_output)
            {
            int abc = I*STEP+J;
            //--------------------------------分子描画(Marvin view)--------------------------------
            int aa=0;
            int Total_1=0, Total_2=0;
            //ラジカル以外
            for(int i=0;i<O;i=i+1){
                for(int j=i+1; j<O; j=j+1){
                    if(B_matrix[i][j]!=0){
                        B_nz[aa].row1=i;
                        B_nz[aa].col1=j;
                        aa=aa+1;
                    }
                }
            }
            Total_1=aa;
            //ラジカル(R-X)
            for(int i=0;i<O;i=i+1){
                for(int j=0; j<B_matrix[i][X]; j=j+1){
                    B_nz[aa].row1=i;
                    aa=aa+1;
                }
            }
            Total_2=aa;


            //以下、Marvin View用の出力ファイル
            fprintf(s,"\n");

            fprintf(s,"  Yoko %d\n",abc);

            fprintf(s,"\n");

            fprintf(s,"%3d%3d  0  0  0  0            999 V2000\n",O+Total_2-Total_1+B_matrix[H][X]+B_matrix[X][X]*2,Total_2+B_matrix[X][X]); //H以外の原子の数，H-R以外の結合の総数

            for (int i=0 ; i<C ;i++){
                fprintf(s,"    0.0       0.0       0.0    C   0  0  0  0  0  0  0  0  0  0  0  0\n");
            }
            for (int i=C ; i<N ; i++){
                fprintf(s,"    0.0       0.0       0.0    N   0  0  0  0  0  0  0  0  0  0  0  0\n");
            }
            for (int i=N ; i<O ; i++){
                fprintf(s,"    0.0       0.0       0.0    O   0  0  0  0  0  0  0  0  0  0  0  0\n");
            }
            for (int i=0 ; i<Total_2-Total_1+B_matrix[H][X]+B_matrix[X][X]*2 ; i++){
                fprintf(s,"    0.0       0.0       0.0    F   0  0  0  0  0  0  0  0  0  0  0  0\n");
            }
            for (int i=0 ; i<Total_1 ; i++){
                fprintf(s,"%3d%3d%3d  0  0  0  0\n",B_nz[i].row1+1,B_nz[i].col1+1,B_matrix[B_nz[i].row1][B_nz[i].col1]);
                /////////①原子の番号１個目、②２個目、③結合次数
            }
            for (int i=Total_1 ; i<Total_2 ; i++){
                fprintf(s,"%3d%3d%3d  0  0  0  0\n",B_nz[i].row1+1,O+1+i-Total_1,1); //R-X
            } 
            for (int i=0 ; i<B_matrix[X][X] ; i++){
                fprintf(s,"%3d%3d%3d  0  0  0  0\n",O+Total_2-Total_1+1+2*i,O+Total_2-Total_1+1+2*i+1,1); //X-X
            } 

            fprintf(s,"M  END\n$$$$\n");

            // abc = I*STEP+J
            }
            #endif
            //--------------------------------分子描画おわり--------------------------------



        } //　J : ステップ
        } // I trial
        free(B_nz); free(R_nz); free(D_energy); free(dH);
    } // omp parallel

    for (int i=0 ; i<PostUV ; i=i+1){
        fprintf(f,"%d %lf %lf %lf\n", i,
                aacid_abundance[i]/(double)TRIAL, sugar_abundance[i]/(double)TRIAL,
                purine_core_abundance[i]/(double)TRIAL);
    }
    for (int i=PostUV ; i<STEP ; i++){
        fprintf(f,"%d %lf %lf %lf\n", i,
                aacid_abundance[i]/(double)TRIAL, sugar_abundance[i]/(double)TRIAL,
                purine_core_abundance[i]/(double)TRIAL);
    }


    fclose(f);
    if(s) fclose(s);
    free(B_energy);
    return 0;
}
                