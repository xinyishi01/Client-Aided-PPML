#include "defines.h"

using namespace std;
using namespace Eigen;

Mat a0(N,D),a1(N,D);
Mat a[2][SIT][5][5], b[2][SIT][5][5], c[2][SIT][5][5];
Mat test_data(testN,D), test_label(testN,1);

int myrandom (int i) { return rand()%i;}

myType randomlong(){
	myType rand1 = abs(rand());
    myType rand2 = abs(rand());
    rand1 = rand1 << (sizeof(int)*8);   
    myType randULL = (rand1 | rand2);   
    return randULL;
}

vector<int> random_perm(){
	vector<int> temp,perm;
	for(int i=0;i<N;i++)
		temp.push_back(i);
	
	for(int i = 0;i<Ep;i++){
		random_shuffle(temp.begin(),temp.end(),myrandom);
		perm.insert(perm.end(),temp.begin(),temp.end());
		
	}
	return perm;
}

void next_batch(Mat& batch,int start, vector<int>& perm, Mat& A){
	for(int i=0;i<B;i++){
		batch.row(i) = A.row(perm[start+i]);
	}
	return;
}

void set_multiplication_triples(vector<int> perm){
	/*
	for(int i=0;i<a0.rows();i++){
		for(int j=0;j<a0.cols();j++){
			a0(i,j) = randomlong();
			a1(i,j) = randomlong();
		}
	}*/

	int lay[4];
	lay[0] = D;
	lay[1] = LAYER1;
	lay[2] = LAYER2;
	lay[3] = 10;

	for (int par = 0; par < 2; ++par)
		for (int i = 0; i < SIT; ++i) {
			for (int j = 1; j <= LL; j++) {
				a[par][i][j][1].resize(B, lay[j - 1]);
				b[par][i][j][1].resize(lay[j - 1], lay[j]);
				c[par][i][j][1].resize(B, lay[j]);
				
				a[par][i][j][2].resize(B, lay[j]);
				b[par][i][j][2].resize(B, lay[j]);
				c[par][i][j][2].resize(B, lay[j]);

				a[par][i][j][3].resize(B, lay[j]);
				b[par][i][j][3].resize(lay[j], lay[j - 1]);
				c[par][i][j][3].resize(B, lay[j - 1]);

				a[par][i][j][4].resize(lay[j - 1], B);
				b[par][i][j][4].resize(B, lay[j]);
				c[par][i][j][4].resize(lay[j - 1], lay[j]);
			} 
		}
	
	for (int par = 0; par < 2; par++)
		for (int i = 0; i < SIT; i++)
			for (int j = 1; j <= LL; j++)
				for (int k = 1; k <= 4; k++) {
					for (int mi = 0; mi < b[par][i][j][k].rows(); mi++)
						for (int mj = 0; mj < b[par][i][j][k].cols(); mj++)
							b[par][i][j][k](mi, mj) = randomlong();
				//	if (j == 1 && k == 1)
				//		continue;
					for (int mi = 0; mi < a[par][i][j][k].rows(); mi++)
						for (int mj = 0; mj < a[par][i][j][k].cols(); mj++)
							a[par][i][j][k](mi, mj) = randomlong();				
				}
	
	for (int par = 0; par < 2; par++)
		for (int i = 0; i < SIT; i++)
			for (int j = 1; j <= LL; j++)
				for (int k = 1; k <= 4; k++) {
				//	if (j == 1 && k == 1)
				//		continue;
				//	cout << par << ' ' << i << ' ' << j << ' ' << k << endl;
					if (k != 2) {
						if (!par) {
							for (int mi = 0; mi < c[par][i][j][k].rows(); mi++)
								for (int mj = 0; mj < c[par][i][j][k].cols(); mj++)
									c[par][i][j][k](mi, mj) = randomlong();
						}
						else {
							c[1][i][j][k] = (a[0][i][j][k] + a[1][i][j][k]) * (b[0][i][j][k] + b[1][i][j][k]) - c[0][i][j][k]; 
						}
					}
					else {
						if (!par) {
							for (int mi = 0; mi < c[par][i][j][k].rows(); mi++)
								for (int mj = 0; mj < c[par][i][j][k].cols(); mj++)
									c[par][i][j][k](mi, mj) = randomlong();
						}
						else {
							for (int mi = 0; mi < c[par][i][j][k].rows(); mi++)
								for (int mj = 0; mj < c[par][i][j][k].cols(); mj++)
									c[par][i][j][k](mi, mj) = (a[0][i][j][k](mi, mj) + a[1][i][j][k](mi, mj)) * (b[0][i][j][k](mi, mj) + b[1][i][j][k](mi, mj)) - c[0][i][j][k](mi, mj); 
						}
					}
				}
/*
	int start_setup = 0;
	for(int i=0;i<SIT;i++){
		Mat c_sum, a_batch0(B,D), a_batch1(B,D), b_sum(D,LAYER1), a_batch(B,D);
		
		next_batch(a_batch0, start_setup, perm,a0);
		next_batch(a_batch1, start_setup, perm,a1);
		a_batch = a_batch0+a_batch1;
		b_sum = b[0][i][1][1] + b[1][i][1][1];
		c_sum = a_batch * b_sum;

		for (int mi = 0; mi < c[0][i][1][1].rows(); mi++)
			for (int mj = 0; mj < c[0][i][1][1].cols(); mj++)
				c[0][i][1][1](mi, mj) = randomlong();
		c[1][i][1][1] = c_sum - c[0][i][1][1]; 
		
		start_setup+=B;
	}*/
/*
	for (int i = 0; i < 2; i++) {
		for (int mi = 0; mi < c[0][i][1][1].rows(); mi++)
			for (int mj = 0; mj < c[0][i][1][1].cols(); mj++) {
				cout << c[0][i][1][1](mi, mj) << ",";
			}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 2; i++) {
		for (int mi = 0; mi < c[0][i][1][1].rows(); mi++)
			for (int mj = 0; mj < c[0][i][1][1].cols(); mj++) {
				cout << c[1][i][1][1](mi, mj) << ",";
			}
		cout << endl;
	}
	cout << endl;*/

	return;
}

void load_train_data(Mat& train_data, Mat& train_label) {
	ifstream infile( "../data/mnist_train.csv" );
	int i=0;
    cout<<"load train data.......\n";
    
	while(infile) {
	//	cout << i << endl;

		string s;
		if (!getline(infile,s))
			break;
		istringstream ss(s);
		int temp;
		char c;
		
		//read label
        ss>>temp;
        ss>>c;
        for (int j = 0; j < 10; j++)
        	train_label(i, j) = 0;
       	train_label(i, temp) = 1;//0~9
       	
        //read data 
        for(int j=0; j<D; j++) {
            ss>>train_data(i,j);
            ss>>c;
        }
		i++;
		if(i>=N) break;
	}

	train_data.conservativeResize(i, D);
	train_label.conservativeResize(i, 10);

	cout<<"n = "<<i<<endl;
	infile.close();
}

int main(int argc, char** argv){
	srand ( unsigned ( time(NULL) ) );
	
	vector<int> perm = random_perm();
	
	cout<<"generating multiplication triples......\n";
	
	int t1=clock();
	
	set_multiplication_triples(perm);
	
	cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s"<<endl;
	
	
	t1=clock();
	
	ofstream F1, F2;
	
	F1.open ("MT0.txt");
	F2.open ("MT1.txt");
	
	for(int i=0;i<perm.size();i++){
		F1<<perm[i]<<",";
		F2<<perm[i]<<",";
		
	}
	F1<<endl;
	F2<<endl;
	/*
	for(int i=0;i<a0.rows();i++){
		for(int j=0;j<a0.cols();j++){
			F1<<a0(i,j)<<",";
			F2<<a1(i,j)<<",";
		}
		F1<<endl;
		F2<<endl;
	}*/

//	for (int par = 0; par < 1; par++)
		for (int j = 1; j <= LL; j++)
			for (int k = 1; k <= 4; k++) {
			//	F1 << j << ' ' << k << endl;
			//	if (j != 1 || k != 1) {
					for (int i = 0; i < SIT; i++)
						for (int mi = 0; mi < a[0][i][j][k].rows(); mi++) {
							for (int mj = 0; mj < a[0][i][j][k].cols(); mj++) {
								F1 << a[0][i][j][k](mi, mj) << ",";
								F2 << a[1][i][j][k](mi, mj) << ",";
							}
							F1 << endl;
							F2 << endl;
						}
			//	}

				for (int i = 0; i < SIT; i++)
					for (int mi = 0; mi < b[0][i][j][k].rows(); mi++) {
						for (int mj = 0; mj < b[0][i][j][k].cols(); mj++) {
							F1 << b[0][i][j][k](mi, mj) << ",";
							F2 << b[1][i][j][k](mi, mj) << ",";
						}
						F1 << endl;
						F2 << endl;
					}

				for (int i = 0; i < SIT; i++)
					for (int mi = 0; mi < c[0][i][j][k].rows(); mi++) {
						for (int mj = 0; mj < c[0][i][j][k].cols(); mj++) {
							F1 << c[0][i][j][k](mi, mj) << ",";
							F2 << c[1][i][j][k](mi, mj) << ",";
						}
						F1 << endl;
						F2 << endl;
					}
			}

	bool t;
    unsigned long rs, r1, r2;
    for (int j = 0; j < IT; j++) {
	    rs = rand();
	    r1 = r2 = 0;
	    for (int i = 63; i >= 0; --i) {
	        if ((rs >> i) & 1) {
	            t = rand() & 1;
	            if (t)
	                r2 += (1 << i);
	            else
	                r1 += (1 << i);
	        }
	    }  
	    F1 << rs << "," << r1 << "," << r2 << ",";
	    F2 << rs << "," << r1 << "," << r2 << ",";
	    F1 << endl;
	  	F2 << endl;
	}
		
	F1.close();
	F2.close();
	
	cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s"<<endl;
	
	t1=clock();
	
	Mat train_data(N,D), train_label(N,10);
	
	load_train_data(train_data, train_label);
	
	cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s"<<endl;
	t1=clock();
	
	cout<<"writing data shares......\n";
	
	Mat train_data0(N,D), train_data1(N,D), train_label0(N,10), train_label1(N,10);
	
	
	for(int i=0;i<train_data.rows();i++){
		for(int j=0;j<train_data.cols();j++){
			train_data0(i,j) = randomlong();
			train_data1(i,j) = train_data(i,j)- train_data0(i,j);
		}
	}
	
	for(int i=0;i<train_label.rows();i++){
		for(int j=0;j<train_label.cols();j++){
			train_label0(i,j) = randomlong();
			train_label1(i,j) = train_label(i,j)- train_label0(i,j);
		}
	}
	
	
	F1.open ("data0.txt");
	F2.open ("data1.txt");
	
	for(int i=0;i<train_data.rows();i++){
		for(int j=0;j<train_data.cols();j++){
			F1<<train_data0(i,j)<<",";
			F2<<train_data1(i,j)<<",";
		}
		F1<<endl;
		F2<<endl;
	}
	
	for(int i=0;i<train_label.rows();i++){
		for(int j=0;j<train_label.cols();j++){
			F1<<train_label0(i,j)<<",";
			F2<<train_label1(i,j)<<",";
		}
		F1<<endl;
		F2<<endl;
	}
	
	F1.close();
	F2.close();
	
	
	/*
	Mat xa= train_data0-a0+train_data1-a1;
	
	F1.open("xa.txt");
	
	for(int i=0;i<train_data.rows();i++){
		for(int j=0;j<train_data.cols();j++){
			F1<<xa(i,j)<<",";
		}
		F1<<endl;
	}
	
	F1.close();
	*/
	cout<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s"<<endl;
}
