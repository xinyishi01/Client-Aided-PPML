#include "utils.hpp"
#include "Layerc.hpp"


class Layer_Client : public Layerc {
public:

	Mat x_mask, W_mask, Go_mask;
	Mat forward_tp, Wg_tp, backward_tp;

	Layer_Client() {}

	Layer_Client(int _D_in, int _D_out) {
		D_in = _D_in;
		D_out = _D_out;
		it = 0;
		x_mask.resize(B, D_in);
		W_mask.resize(D_in, D_out);
		Go_mask.resize(B, D_out);
		forward_tp.resize(B, D_out);
		Wg_tp.resize(D_in, D_out);
		backward_tp.resize(B, D_in);
		forward_z1.resize(B, D_out);
		forward_z2.resize(B, D_out);
		
		for (int i = L; i >= 0; i--)
			b_mask[i].resize(B, D_out);
		
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
    
    Mat *nwX = new Mat[rM];
    Mat *nwW = new Mat[rM];
    Mat *nwG = new Mat[rM];
    
    for (int i = 0; i < rM; ++i) {
        nwX[i].resize(B, D_in); 
        nwW[i].resize(D_in, D_out);
        nwG[i].resize(B, D_out);
        add_mat(prf0, prf1, nwX[i]); 
        add_mat(prf0, prf1, nwW[i]); 
        add_mat(prf0, prf1, nwG[i]);
    }
    
	Mat temp0(rM, B * D_out), temp1(rM, B * D_in), temp2(rM, D_in * D_out);
    
    for (int i = 0; i < rM; ++i) {
    	temp0.row(i) = flatten_mat(nwX[i] * nwW[i]);
    	temp1.row(i) = flatten_mat(nwG[i] * nwW[i].transpose());
    	temp2.row(i) = flatten_mat(nwX[i].transpose() * nwG[i]);
    }
    
    distribute_mat_same(temp0, prf0, prf1, io_alice, io_bob);
    distribute_mat_same(temp1, prf0, prf1, io_alice, io_bob);
    distribute_mat_same(temp2, prf0, prf1, io_alice, io_bob);
    
    delete prf0;
    delete prf1;
    
    delete[] nwX;
    delete[] nwW;
    delete[] nwG;

}

/*--------------------------Malicious--------------------------*/	


	void pre(PRF *_prf_alice, PRF *_prf_bob, PRF *_prf_alice_same, PRF *_prf_bob_same) {
		/*--------------------------Malicious--------------------------*/
		
		
		if(oNum == 0) {
			//cout << "L2/L3 is generating triples...\n";
			check_malicious();
			oNum = 0;
		}
		
		if((++oNum) == rN) oNum = 0;
		/*--------------------------Malicious--------------------------*/
		
		prf_alice = _prf_alice;
		prf_bob = _prf_bob;
		prf_alice_same = _prf_alice_same;
		prf_bob_same = _prf_bob_same;
		
		add_mat2(prf_alice, prf_bob, forward_z1);
		if(id != "L3") add_mat2(prf_alice, prf_bob, forward_z2);
		
		if (id == "L3") {
			for (int i = L; i >= 0; i--)
				add_mat2(prf_alice, prf_bob, b_mask[i]);
		}
		
		
	}
	
};
