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
	}
	
	void pre(PRF *_prf_alice, PRF *_prf_bob, PRF *_prf_alice_same, PRF *_prf_bob_same) {
		prf_alice = _prf_alice;
		prf_bob = _prf_bob;
		prf_alice_same = _prf_alice_same;
		prf_bob_same = _prf_bob_same;
		
		add_mat2(prf_alice, prf_bob, x_mask);
		add_mat(prf_alice_same, prf_bob_same, W_mask);
		add_mat2(prf_alice, prf_bob, Go_mask);
		add_mat2(prf_alice, prf_bob, forward_z1);
		if(id != "L3") add_mat2(prf_alice, prf_bob, forward_z2);

		if (id == "L3") {
			for (int i = L; i >= 0; i--)
				add_mat2(prf_alice, prf_bob, b_mask[i]);
		}
		
		forward_tp = x_mask * W_mask;
		backward_tp = Go_mask * (W_mask.transpose());
		
		Mat *temp = new Mat(B, D_in * D_out);
		for (int k = 0; k < B; ++k) {
			Wg_tp = (x_mask.row(k)).transpose() * (Go_mask.row(k));
			for (int i = 0, cnt = 0; i < D_in; ++i) for (int j = 0; j < D_out; ++j, ++cnt) {
				(*temp)(k, cnt) = Wg_tp(i, j);
			}
		}
		
		distribute_mat(forward_tp, prf_alice, prf_bob, io_alice, io_bob);
		distribute_mat(backward_tp, prf_alice, prf_bob, io_alice, io_bob);
		distribute_mat(*temp, prf_alice, prf_bob, io_alice, io_bob);
		
		delete temp;
		
	}
	
};
