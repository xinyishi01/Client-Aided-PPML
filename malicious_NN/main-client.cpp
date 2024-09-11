#include "utils.hpp"
#include "NN_Client.hpp"

NetIO *io_server;
NetIO *io_client_alice;
NetIO *io_client_bob;
NetIO *io_client;

void setup(int party) {
    if (party == ALICE || party == EVE) {
        io_client_alice = new NetIO(party==EVE ? nullptr : "127.0.0.1", 1027);
        io_client_alice->set_nodelay();
    }
    if (party == BOB || party == EVE) {
        io_client_bob = new NetIO(party==EVE ? nullptr : "127.0.0.1", 1028);
        io_client_bob->set_nodelay();
    }
	
	if(party == ALICE) io_client = io_client_alice;
	if(party == BOB) io_client = io_client_bob;
}

void load_train_data(Mat& train_data, Mat& train_label) {
	ifstream infile( "../data/mnist_train.csv" );
	int i=0;
    cout<<"load train data.......\n";
    
    train_data.resize(N, D);
    train_label.resize(N, K);
    for (int i = 0; i < N; ++i) for (int j = 0; j < K; ++j) train_label(i, j) = 0;
    
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
        train_label(i, temp) = 1 << L;
        
        //read data 
        for(int j=0; j<D; j++) {
            ss>>train_data(i,j);
            train_data(i, j) <<= L;
            ss>>c;
        }
		i++;
		if(i>=N) break;
	}
	
	if(i < N) assert(0);


	cout<<"n = "<<i<<endl;
	infile.close();
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


int main(int argc, char** argv) {
	int port, party;
	parse_party_and_port(argv, &party, &port);
	setup(party);
	
	Mat train_data(N,D), train_label(N,K);
	load_train_data(train_data, train_label);
	
	vector<int> perm = random_perm();
	io_client_alice->send_data(&perm[0], sizeof(int) * perm.size());
	io_client_alice->flush();
	io_client_bob->send_data(&perm[0], sizeof(int) * perm.size());
	io_client_bob->flush();
	
	NN_Client nn_client(io_client_alice, io_client_bob);
	
	nn_client.online(train_data, train_label, perm);
	return 0;
}

