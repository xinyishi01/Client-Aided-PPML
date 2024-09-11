#include "utils.hpp"
#include "Layerc.hpp"


class Layer1_Client : public Layerc {
public:

	Mat z, tau, backward_z, backward_tau;
	
	Layer1_Client() {}

	Layer1_Client(int _D_in, int _D_out) {
		D_in = _D_in;
		D_out = _D_out;
		it = 0;
		z.resize(D_in, D_out);
		tau.resize(B, D_out);
		backward_z.resize(B, D_out);
		backward_tau.resize(D_in, D_out);
		forward_z1.resize(B, D_out);
		forward_z2.resize(B, D_out);
		
		/*--------------------------Malicious--------------------------*/
		oNum = 0;
		dNum = 0;
		/*--------------------------Malicious--------------------------*/
	}

/*--------------------------Malicious--------------------------*/	

void check_malicious() {
    PRF *prf0 = new PRF, *prf1 = new PRF;
    receive_keys(prf0, io_alice, 1);
    receive_keys(prf1, io_bob, 1);
    
    Mat *nwX = new Mat[cM];
    Mat *nwZ = new Mat[cM];
    Mat *nwbZ = new Mat[cM];
    
    for (int i = 0; i < cM; ++i) {
        nwX[i].resize(B, D_in); 
        nwZ[i].resize(D_in, D_out);
        nwbZ[i].resize(B, D_out);
        add_mat(prf0, prf1, nwX[i]); 
        add_mat(prf0, prf1, nwZ[i]); 
        add_mat(prf0, prf1, nwbZ[i]);
    }
    
    Mat temp0(cM, B * D_out), temp1(cM, D_in * D_out);
    
    for (int i = 0; i < cM; ++i) {
    	temp0.row(i) = flatten_mat(nwX[i] * nwZ[i]);
    	temp1.row(i) = flatten_mat(nwX[i].transpose() * nwbZ[i]);
    }
    
    distribute_mat_same(temp0, prf0, prf1, io_alice, io_bob);
    distribute_mat_same(temp1, prf0, prf1, io_alice, io_bob);
    
    delete prf0;
    delete prf1;
    
    delete[] nwX;
    delete[] nwZ;
    delete[] nwbZ;

}

/*--------------------------Malicious--------------------------*/


	void pre(Mat in_x, PRF *_prf_alice, PRF *_prf_bob, PRF *_prf_alice_same, PRF *_prf_bob_same) {
		prf_alice = _prf_alice;
		prf_bob = _prf_bob;
		prf_alice_same = _prf_alice_same;
		prf_bob_same = _prf_bob_same;
	
		add_mat(prf_alice_same, prf_bob_same, z);
		add_mat2(prf_alice, prf_bob, forward_z1);
		add_mat2(prf_alice, prf_bob, forward_z2);
		add_mat2(prf_alice, prf_bob, backward_z);
		tau = in_x * z;
		Mat *temp = new Mat(B, D_in * D_out);
		
		for (int k = 0; k < B; ++k) {
			backward_tau = (in_x.row(k)).transpose() * (backward_z.row(k));
			for (int i = 0, cnt = 0; i < D_in; ++i) for (int j = 0; j < D_out; ++j, ++cnt) {
				(*temp)(k, cnt) = backward_tau(i, j);
			}
		}
		
		distribute_mat(tau, prf_alice, prf_bob, io_alice, io_bob);
		distribute_mat(*temp, prf_alice, prf_bob, io_alice, io_bob);
		
		delete temp;
		
		/*--------------------------Malicious--------------------------*/
		
		++oNum;
		
		if(oNum == cN) {
			//cout << "L1 is checking malicious...\n";
			check_malicious();
			oNum = 0;
		}
		/*--------------------------Malicious--------------------------*/
		
	}

};
