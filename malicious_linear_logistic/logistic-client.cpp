#include "utils.hpp"

typedef uint8_t* vec_t;

NetIO * io_server;
NetIO * io_client_alice;
NetIO * io_client_bob;
NetIO * io_client;

PRF *prf_alice, *prf_bob, *prf_alice_same, *prf_bob_same;

/*-------------------Malicious----------------*/

int cNum0(0);
Mat cX0[cM + 10], cY0[cM + 10];
Mat cX1[cM + 10], cY1[cM + 10];

/*-------------------Malicious----------------*/

/*-------------------Malicious-act----------------*/

int dNum(0);

/*-------------------Malicious-act----------------*/

void load_train_data(Mat& train_data, Mat& train_label){
    ifstream infile( "../data/mnist_train.csv" );
        int count1=0, count2=0;
        int i=0;
        while(infile) {

            string s;
            if (!getline(infile,s))
                break;
            istringstream ss(s);
            int temp;
            char c;


            //read label
            ss>>temp;
            ss>>c;
            if(temp == 0 && count1<N/2) {
                train_label(i) = 0;
                count1++;

                //read data (last entry 1)
                for(int j=0; j<D-1; j++) {
                    ss>>train_data(i,j);
                    ss>>c;
                }

                train_data(i,D-1) = 1;
                i++;
            }


            if(temp != 0 && count2<N/2) {
                train_label(i) = 1;
                count2++;

                //read data (last entry 1)
                for(int j=0; j<D-1; j++) {
                    ss>>train_data(i,j);
                    ss>>c;
                }

                train_data(i,D-1) = 1;
                i++;
            }


            if(i>=N)
                break;
        }

        //train_data.conservativeResize(i, D);
        //train_label.conservativeResize(i,1);
    infile.close();
}

static bool comp(const vec_t &repr1, const vec_t &repr2) {
	for (int i = 0; i < 16; i++) {
		if (repr1[i] < repr2[i]) return 1;
		if (repr1[i] > repr2[i]) return 0;
	}
	return 0;
}

bool check_reprs_equal(vec_t repr1, vec_t repr2, int length) {
	for (int i = 0; i < length; i++) {
		if (repr1[i] != repr2[i]) return false;
	}
	return true;
}

bool compute_vec_intersection(vec_t* vec1, vec_t* vec2) {
	unsigned int count = 0;
	int j = 0;
	for (int i = 0; i < P; i++){
		while (j < P && comp(vec2[j], vec1[i])) ++j;
		if (j >= P) break;
		if (check_reprs_equal(vec1[i], vec2[j], sizeof(unsigned long) * 2)) {
			count++;
			//cerr << i << ' ' << j << ' ' << vec1[i][0] << endl;
		}
	}
	assert(count == 0 || count == 1);
	return (bool)count;
}

void receive_encodings(Mat &counts1, Mat &counts2, NetIO* alice_io, NetIO* bob_io) {
	//vector<uint8_t> enc_t0(2 * B * P * 16, 0), enc_t1(2 * B * P * 16, 0);
	uint8_t enc_t0[2 * B * P * 16], enc_t1[2 * B * P * 16];
	alice_io->recv_data(&enc_t0[0], sizeof(uint8_t) * 2 * B * P * 16);
    bob_io->recv_data(&enc_t1[0], sizeof(uint8_t) * 2 * B * P * 16);
    
    int cnt0 = 0, cnt1 = 0;
    for (int i = 0; i < 2; i++)
    	for (int j = 0; j < B; j++){
    		vec_t alice_enc[P], bob_enc[P];

    		for (int k = 0; k < P; k++){
    			alice_enc[k] = new uint8_t[sizeof(unsigned long) * 2];
    			for (int l = 0; l < 16; l++)
    				alice_enc[k][l] = enc_t0[cnt0++];
    		}
    		for (int k = 0; k < P; k++){
    			bob_enc[k] = new uint8_t[sizeof(unsigned long) * 2];
    			for (int l = 0; l < 16; l++)
    				bob_enc[k][l] = enc_t1[cnt1++];
    		}
    		
    		if (i == 0)
    			counts1(j, 0) = compute_vec_intersection(alice_enc, bob_enc);
    		else
    			counts2(j, 0) = compute_vec_intersection(alice_enc, bob_enc);
    		
    		for (int k = 0; k < P; k++) {
				delete[] alice_enc[k];
			}
			for (int k = 0; k < P; k++) {
				delete[] bob_enc[k];
			}
    	}
}


void activation_function_client(Mat &z1, Mat &z2, int fl) {
	Mat counts1(B, 1), counts2(B, 1);
	Mat tau1(B, 1), tau2(B, 1);
	
	receive_encodings(counts1, counts2, io_client_alice, io_client_bob);
//	cerr << "aft rec\n";

	tau1 = counts1.cwiseProduct(z1);
	tau2 = counts2.cwiseProduct(z2);

//	cerr << "bef\n";
	
	unsigned long int idx_alice = prf_alice_same->index;
	unsigned long int idx_bob = prf_bob_same->index;
	distribute_mat(counts1, prf_alice, prf_bob, io_client_alice, io_client_bob);
	distribute_mat(counts2, prf_alice, prf_bob, io_client_alice, io_client_bob);
	distribute_mat(tau1, prf_alice, prf_bob, io_client_alice, io_client_bob);
	distribute_mat(tau2, prf_alice, prf_bob, io_client_alice, io_client_bob);
	if (fl) {
        prf_alice_same->index = idx_alice;
        prf_bob_same->index = idx_bob;
	}
}



/*-------------------Malicious----------------*/
void check_malicious() {
//	return;
    PRF *prf0 = new PRF, *prf1 = new PRF;
    receive_keys(prf0, io_client_alice, 1);
    receive_keys(prf1, io_client_bob, 1);
    
    for (int i = 0; i < cM; ++i) {
    	cX0[i].resize(B, D); cY0[i].resize(D, 1);
    	add_mat(prf0, prf1, cX0[i]); add_mat(prf0, prf1, cY0[i]);
    	/*cX1[i].resize(D, 1);*/ cY1[i].resize(1, B);
    	/*add_mat(prf0, prf1, cX1[i]);*/ add_mat(prf0, prf1, cY1[i]);
    }
    
    Mat temp0(B, cM), temp1(cM, D);

    Mat aa(1, D);
    temp1.row(0) = aa;
    
    for (int i = 0; i < cM; ++i) {
    	temp0.col(i) = cX0[i] * cY0[i];
    	temp1.row(i) = sp_product(cX0[i], cY1[i]);
    }
    
    distribute_mat_same(temp0, prf0, prf1, io_client_alice, io_client_bob);
    distribute_mat_same(temp1, prf0, prf1, io_client_alice, io_client_bob);
    
    delete prf0;
    delete prf1;

}

/*-------------------Malicious----------------*/

/*-------------------Malicious-act----------------*/

void receive_check_encodings(Mat &counts, NetIO* alice_io, NetIO* bob_io, int dM) {

	uint8_t* enc_t0 = new uint8_t[dM * P * 16];
	uint8_t* enc_t1 = new uint8_t[dM * P * 16];
	alice_io->recv_data(&enc_t0[0], sizeof(uint8_t) * dM * P * 16);
    bob_io->recv_data(&enc_t1[0], sizeof(uint8_t) * dM * P * 16);

    int cnt0 = 0, cnt1 = 0;
    for (int i = 0; i < 1; i++)
    	for (int j = 0; j < dM; j++){
    		vec_t alice_enc[P], bob_enc[P];

    		for (int k = 0; k < P; k++){
    			alice_enc[k] = new uint8_t[sizeof(unsigned long) * 2];
    			for (int l = 0; l < 16; l++)
    				alice_enc[k][l] = enc_t0[cnt0++];
    		}
    		for (int k = 0; k < P; k++){
    			bob_enc[k] = new uint8_t[sizeof(unsigned long) * 2];
    			for (int l = 0; l < 16; l++)
    				bob_enc[k][l] = enc_t1[cnt1++];
    		}
    		
    		counts(j, 0) = compute_vec_intersection(alice_enc, bob_enc);
    		
    		for (int k = 0; k < P; k++) { delete[] alice_enc[k]; } for (int k = 0; k
< P; k++) { delete[] bob_enc[k]; } }

    delete []enc_t0;
    delete []enc_t1;
}

void check_malicious_activation_function(int dN) {
	int dM = dN * dk + dk + 1;
	Mat er(dM, 1);
	add_mat(prf_alice_same, prf_bob_same, er);

	Mat counts(dM, 1), tau(dM, 1);

	receive_check_encodings(counts, io_client_alice, io_client_bob, dM);
	tau = counts.cwiseProduct(er);
	
	distribute_mat_same(counts, prf_alice_same, prf_bob_same, io_client_alice, io_client_bob);
	distribute_mat_same(tau, prf_alice_same, prf_bob_same, io_client_alice, io_client_bob);
}

/*-------------------Malicious-act----------------*/

int main(int argc, char** argv){	
	//srand((unsigned)time(NULL));

    int port, party;
    parse_party_and_port(argv, &party, &port);
    setup(party, io_server, io_client_alice, io_client_bob, io_client);

    cout << "reading data......\n";

    Mat train_data(N,D), train_label(N,1);

     load_train_data(train_data, train_label);

    vector<int> perm = random_perm();
    
    io_client_alice->send_data(&perm[0], sizeof(int) * perm.size());
    io_client_bob->send_data(&perm[0], sizeof(int) * perm.size());
   
    int hahaha1;
    io_client_alice->recv_data(&hahaha1, sizeof(int));
	io_client_bob->recv_data(&hahaha1, sizeof(int));
	
//	cerr << "11\n";
	
   
	int start = 0;

    
    Mat x_batch(B, D), y_batch(B, 1), tau0_batch(B, 1), tau1_batch(B, D);
    Mat z0(1, D), z1(B, 1);
    prf_alice = new PRF[B];
    prf_bob = new PRF[B];
    prf_alice_same = new PRF;
    prf_bob_same = new PRF;
    Mat z0_act_batch(B, 1), z1_act_batch(B, 1);
    
    cout<<"start training..."<<endl;
    
    PRF *prf_alice_temp = new PRF[N];
	PRF *prf_bob_temp = new PRF[N];
	
	receive_keys(prf_alice_temp, io_client_alice, N);
	receive_keys(prf_bob_temp, io_client_bob, N);
	
	distribute_mat(train_data, prf_alice_temp, prf_bob_temp, io_client_alice, io_client_bob);
	distribute_mat(train_label, prf_alice_temp, prf_bob_temp, io_client_alice, io_client_bob);
	
	delete [] prf_alice_temp;
	delete [] prf_bob_temp;
    
	
	//start training
	for (int i = 0; i < IT; ++i, start += B) {
        
    	next_batch(x_batch,start,perm,train_data);
        next_batch(y_batch,start,perm,train_label);
        
        //cout << "receive keys...\n";
        receive_keys(prf_alice, io_client_alice, B);
        receive_keys(prf_bob, io_client_bob, B);
        
        receive_keys(prf_alice_same, io_client_alice, 1);
        receive_keys(prf_bob_same, io_client_bob, 1);
        
    //    cout << "add mat...\n";
        
        add_mat(prf_alice_same, prf_bob_same, z0);
        add_mat2(prf_alice, prf_bob, z1);
        
    //    cout << "calculation...\n";
        
        for (int j = 0; j < B; ++j) {
        	tau0_batch(j, 0) = ((x_batch.row(j)).cwiseProduct(z0)).sum();
        	//tau1_batch.row(j) += x_batch.row(j) * z1(j, 0);
        	for (int k = 0; k < D; ++k) tau1_batch(j, k) = x_batch(j, k) * z1(j, 0);
        }
        
    //    cout << "distribute mat...\n";
		distribute_mat(tau0_batch, prf_alice, prf_bob, io_client_alice, io_client_bob);
		distribute_mat(tau1_batch, prf_alice, prf_bob, io_client_alice, io_client_bob);
		
		if(DEBUG_PRF) {
			send_mat(tau1_batch, io_client_alice);
		}
		
		add_mat(prf_alice_same, prf_bob_same, z0_act_batch);
		add_mat(prf_alice_same, prf_bob_same, z1_act_batch);
		
        activation_function_client(z0_act_batch, z1_act_batch, 0);
        //cout << "activation function done ...\n";
        
        /*-------------------Malicious----------------*/
		cNum0 += 1;
		
		if (cNum0 == cN) {
            check_malicious();
            cNum0 = 0;
        }

        /*-------------------Malicious----------------*/

        /*-------------------Malicious-act----------------*/
        dNum += 2 * B;

        if (dNum >= dNN) {
        //	cerr << "check_malicious_activation_function\n";
        	check_malicious_activation_function(dNum);
        	dNum = 0;
        }
        /*-------------------Malicious-act----------------*/
    }

	delete[] prf_alice;
	delete[] prf_bob;
	delete prf_alice_same;
	delete prf_bob_same;
	
	if (party != EVE) {
		delete io_server;
	}
	if (party != BOB) {
		delete io_client_alice;
	}
	if (party != ALICE) {
		delete io_client_bob;
	}
}
