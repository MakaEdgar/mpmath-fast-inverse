#include <iostream>
#include <chrono>
using namespace std::chrono;


#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <Eigen/Dense>
const int MP_PRECISION = 50;
typedef boost::multiprecision::
number<boost::multiprecision::backends::cpp_dec_float<MP_PRECISION>> doubleMP;
typedef Eigen::Matrix<doubleMP, Eigen::Dynamic, Eigen::Dynamic> MatrixMP;


MatrixMP& _get_matrix2(void* ptr_mpmat) {
	assert(ptr_mpmat && "matrix is not initialized");
	return *((MatrixMP*)ptr_mpmat);
}




#include "mpinv.h"


int main(int, char**)
{
	using namespace Eigen;
	std::cout << std::setprecision(std::numeric_limits<doubleMP>::digits10);

	const char* mpmat_file = "8x8.mpmat";
	const char* mpmat_inv_file = "8x8_inv.mpmat";

	mpmat m(mpmat_file);
	m.calc_invert_with_logdet(true);
	m.save_mpmat_inv(mpmat_inv_file);
	
	mpmat m_inv(mpmat_inv_file);
	m_inv.calc_invert_with_logdet(true);
	m_inv.save_mpmat_inv(mpmat_file);

	mpmat m_inv_inv(mpmat_file);
	m_inv_inv.calc_logdet(true);
	
	std::cout << "Matrix logdet 1 is " << m.get_logdet() << std::endl;
	std::cout << "Matrix logdet 2 is " << m_inv.get_logdet() << std::endl;
	std::cout << "Matrix logdet 3 is " << m_inv_inv.get_logdet() << std::endl;


	const char *data = "8 8 \
		12.731464353632530484679290238111378727041055499879 12.724763332504386521579960479434732462487578072276 12.716995152109991605717120986502562657298385633336 12.707924281955546626144117509686738947706958872409 12.697842392248626924436482642775755228448156590264 12.735365901619125238834264469986692285286168795112 12.690799787932745982797919103506234361878253448296  12.71024547547631388762261591474889442294810471193 \
		12.724763332504386521579960479434732462487578072276 12.722566879950415617373972674603816484795146716015 12.719061765406726394768140172281776209073194415966 12.714131617841576857999057500156474336739801622439 12.707924281955546626144117509686738947706958872409 12.725480920235890843928194113243942115368403688929 12.703247117957455022149542269380260818739083824587 12.715461466279509144729164262162925823713365395349 \
		12.716995152109991605717120986502562657298385633336 12.719061765406726394768140172281776209073194415966 12.719722299556828895376904807259237934328295623536 12.718974601408975078186381594121954935455084524072 12.716821160862515866584743392947481482176625725026 12.714779093715384769787662508719777622253200683332 12.714696825311724751121432638251850073497568603657 12.719289902295858727018862679167554139720975263445 \
		12.707924281955546626144117509686738947706958872409 12.714131617841576857999057500156474336739801622439 12.718974601408975078186381594121954935455084524072 12.722566879950415617373972674603816484795146716015 12.724763332504386521579960479434732462487578072276 12.702939800490084788917455538521243622491819035218 12.725456510264436294787793984432216425468488979363 12.721816052282085859280215222598945582929412660611 \
		12.697842392248626924436482642775755228448156590264 12.707924281955546626144117509686738947706958872409 12.716821160862515866584743392947481482176625725026 12.724763332504386521579960479434732462487578072276 12.731464353632530484679290238111378727041055499879 12.690349247765232098184510919396771154435102854936 12.735153326167135464378891770425089022292084456839 12.722930573571887042974451382001027820952986285226 \
		12.735365901619125238834264469986692285286168795112 12.725480920235890843928194113243942115368403688929 12.714779093715384769787662508719777622253200683332 12.702939800490084788917455538521243622491819035218 12.690349247765232098184510919396771154435102854936 12.741567273757554465364793754654078802232590331761 12.681819869156522719450457790806545575687073736917 12.705916029786194315835917047441281917338955155529 \
		12.690799787932745982797919103506234361878253448296 12.703247117957455022149542269380260818739083824587 12.714696825311724751121432638251850073497568603657 12.725456510264436294787793984432216425468488979363 12.735153326167135464378891770425089022292084456839 12.681819869156522719450457790806545575687073736917 12.740876267442460549224911876700972938970552147247 12.722915815133745633192263381794760505869916848545 \
		12.71024547547631388762261591474889442294810471193 12.715461466279509144729164262162925823713365395349 12.719289902295858727018862679167554139720975263445 12.721816052282085859280215222598945582929412660611 12.722930573571887042974451382001027820952986285226 12.705916029786194315835917047441281917338955155529 12.722915815133745633192263381794760505869916848545 12.721328362066752123580300159606306018495780724566";

	m.get_data(data);




	//mpmat m(in_mpmat_file);
	//m.save_mpmat("qqq.mpmat");
	//m.calc_invert();
	//m.save_mpmat_inv("qqq_inv.mpmat");


	//mpmat m1(in_mpmat_file);
	//mpmat m2(in_mpmat_file);
	//mpmat m3(in_mpmat_file);

	//auto start = high_resolution_clock::now();
	//m1.logdet = 0;
	//m1.calc_logdet1(true);

	//auto middle1 = high_resolution_clock::now();
	//m2.logdet = 0;
	//m2.calc_logdet2(true);

	//auto middle2 = high_resolution_clock::now();
	//m3.logdet = 0;
	//m3.calc_logdet2(true);

	//auto stop = high_resolution_clock::now();

	//std::cout << "logdet1 " << duration_cast<milliseconds>(middle1 - start).count() << std::endl;
	//std::cout << "logdet2 " << duration_cast<milliseconds>(middle2 - middle1).count() << std::endl;
	//std::cout << "logdet3 " << duration_cast<milliseconds>(stop - middle2).count() << std::endl;

	//std::cout << std::endl;

	//std::cout << "Matrix logdet 1 is " << m1.logdet << std::endl;
	//std::cout << "Matrix logdet 2 is " << m2.logdet << std::endl;
	//std::cout << "Matrix logdet 3 is " << m3.logdet << std::endl;


	//mpmat m_inv = mpmat();
	//m_inv.load_mpmat("qqq_inv.mpmat");
	//m_inv.calc_logdet2(true);
	//std::cout << "Matrix inv logdet2 is " << m_inv.logdet << std::endl;



	////MatrixMP C = MatrixMP(*((MatrixMP*)m.ptr_mpmat));
	//MatrixMP C = _get_matrix2(m.ptr_mpmat);

	//doubleMP d = C.determinant();


	//m.calc_invert_with_logdet(true);




	


	return 0;
}










// old main
// int main(int, char**) {
//int N = 20;
//MatrixMP A = MatrixMP::Random(N, N);
//MatrixMP V(A);
//V.transposeInPlace();
//A += V;

//std::cout << "Here is the matrix A:\n" << A << std::endl;
//std::cout << std::endl;
//std::cout << "The inverse of A is:\n" << A.inverse() << std::endl;

//std::cout << "Start inversing matrix!" << std::endl;
//MatrixMP B = A.inverse();
//std::cout << "End inversing matrix!" << std::endl;
//std::cout << std::endl;

//std::ofstream ofs("inverse.txt");
//ofs << std::setprecision(std::numeric_limits<doubleMP>::digits10);
//ofs << B << std::endl;



//doubleMP detA = A.determinant();
//doubleMP detB = B.determinant();


//std::cout << "Matrix determinant is " << detA << std::endl;
//std::cout << "Inv Matrix determinant is " << detB << std::endl;
//std::cout << "Mat * InvMat det product is " << detA * detB << std::endl;
//std::cout << std::endl;

//double detA_ = detA.convert_to<double>();
//std::cout << "detA to double is " << detA_ << std::endl;
//std::cout << "detA double log is " << std::log(detA_) << std::endl;
//std::cout << std::endl;

//std::cout << "Alternative logdet: " << logdet(A, false) << std::endl;
//std::cout << "Alternative logdet2: " << logdet(A, true) << std::endl;
//}