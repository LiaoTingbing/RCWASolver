#include "Device.h"
#include "userFunc.h"



Device::Device()
{

	string filename[] = {
"Izpos",
"ku",
"kv",
"lambda",
"theta",
"phi",
"x",
"y",
"z",
"IndexZ"
	};
	for (auto& p : filename)
	{
		p = "input/" + p + ".txt";
	}

	Izpos.load(filename[0]);
	loadTXT(ku, filename[1]);
	loadTXT(kv, filename[2]);
	loadTXT(lambda, filename[3]);
	loadTXT(theta, filename[4]);
	loadTXT(phi, filename[5]);
	x.load(filename[6]);
	y.load(filename[7]);
	z.load(filename[8]);


	layersNum = Izpos.size();

	string s;
	IndexZ = new cx_mat[layersNum];
	for (int i = 0; i < layersNum;i++)
	{
		s = "input/IndexZ" + to_string(i + 1) + ".txt";
		//cout << s << endl	;
		IndexZ[i].load(s);
		//cout << IndexZ[i].n_rows << IndexZ[i].n_cols << endl;
	}



}

Device::~Device()
{
	delete[]   IndexZ;
}

void Device::RCWA()
{
	//cout << "\t ÖÐÐÄÎ»ÖÃ" << Izpos.st() << endl;
	//cout << "\t ku" << ku << endl;
	//cout << "\t kv" << kv << endl;
	//cout << "\t lambda" << lambda << endl;
	//cout << "\t theta" << theta << endl;
	//cout << "\t phi" << phi << endl;
	//cout << "\t x" << x.st() << endl;
	//cout << "\t y" << y.st() << endl;
	//cout << "\t z" << z.st() << endl;
	//
	Nx = x.size();
	Ny = y.size();


	Ncx = x.size() - 1;
	Ncy = y.size() - 1;

	Nc = Ncx * Ncy;

	Lx = x(Nx - 1) - x(0);
	Ly = y(Ny - 1) - y(0);
	thetaArc = theta * pi / 180.0;
	phiArc = phi * pi / 180.0;
	vec incident_direction = {
		sin(thetaArc) * cos(phiArc),
		sin(thetaArc) * sin(phiArc),
		cos(thetaArc)
	};
	double k0 = 2 * pi / lambda;


	cx_double ER_inc = 1.4 * 1.4;
	cx_double UR_inc = 1;

	cx_vec k_inc = sqrt(ER_inc) * incident_direction;

	cx_double kx_inc = k_inc(0);
	cx_double ky_inc = k_inc(1);
	cx_double kz_inc = k_inc(2);

	double Tx = 2.0 * pi / Lx / k0;
	double Ty = 2.0 * pi / Ly / k0;

	ku = ku;
	kv = kv;

	vec m = linspace(-ku, ku, 2 * ku + 1);
	vec n = linspace(-kv, kv, 2 * kv + 1);

	size_t M = m.size();
	size_t N = n.size();


	cx_mat kx_mn(M, N, fill::zeros);
	for (int i = 0; i < kx_mn.n_cols;i++)
	{
		kx_mn.col(i) = m * Tx + kx_inc;
	}
	cx_mat ky_mn(M, N, fill::zeros);
	for (int i = 0; i < ky_mn.n_rows;i++)
	{
		ky_mn.row(i) = n.st() * Ty + ky_inc;
	}

	//ky_mn.print();

	cx_mat Kx = matrix_to_vector_row(kx_mn);
	cx_mat Ky = matrix_to_vector_row(ky_mn);

	size_t Nh = M * N;

	cx_mat I(Nh, Nh, fill::eye);
	cx_mat II(2 * Nh, 2 * Nh, fill::eye);
	cx_mat Z(Nh, Nh, fill::zeros);
	cx_mat ZZ(Nh * 2, Nh * 2, fill::zeros);


	G.S11 = ZZ;
	G.S12 = II;
	G.S21 = II;
	G.S22 = ZZ;


	WVLmatrix W0V0 = Gap(Kx, Ky, I, Z);

	cx_mat W0 = W0V0.W;
	cx_mat V0 = W0V0.V;

	cx_double  ER_ref = ER_inc;
	cx_double  UR_ref = UR_inc;

	WVLmatrix WVref = solution_in_Homogeneous_Layers(Kx, Ky, ER_ref, UR_ref);
	cx_mat Wref = WVref.W;
	cx_mat Vref = WVref.V;

	cx_mat Aref = inv(W0) * Wref + inv(V0) * Vref;
	cx_mat Bref = inv(W0) * Wref - inv(V0) * Vref;

	Smatrix SR;
	SR.S11 = -inv(Aref * Bref);
	SR.S12 = 2.0 * II * inv(Aref);
	SR.S21 = 0.5 * (Aref - Bref * inv(Aref) * Bref);
	SR.S22 = Bref * inv(Aref);

	G = SconnectRight(G, SR);

	// Z cell

	vec ZL = diff(z);
	ZL.print();
	cout << ZL.size();

	double di;
	size_t NumIndex;
	cx_mat ERi, unERi;

	for (size_t Layer = 0; Layer < ZL.size();Layer++)
	{
		di = ZL(Layer);
		ERi = pow(IndexZ[Layer], 2.0);

		unERi = unique(ERi);
		//cout << unERi.size() << endl;;
		if (unERi.size() == 1)
		{

		}

	}







}


