#ifndef MASTER_H
#define MASTER_H 

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <iomanip>
#include <SymEigsSolver.h> 
#include <Eigen/Eigenvalues>
#include <sstream>
#include <fstream>


using namespace Eigen;
using namespace std;
using namespace Spectra;

class master
{
private:
        MatrixXcd _rho;
        MatrixXd a, adag;
        MatrixXd sigmamin, sigmaplus, sigmaz;
        complex<double> i;
public:
        static int Max;

        const MatrixXcd& rho()const{return _rho;};

        master(const double& deltar, const double& deltaq, const double& j,
                const double& f, const double& kappa, const double& gamma, ofstream& outfile):
        _rho(MatrixXcd::Zero(2*(Max+1), 2*(Max+1))),
        i(0,1)
        {
                Matrix2d tempsigmaz, tempsigmamin, tempsigmaplus, tempsigmaeye;
                tempsigmaz<<1,0,0,-1;tempsigmamin<<0,0,1,0;tempsigmaplus<<0,1,0,0;tempsigmaeye<<1,0,0,1;
                MatrixXd tempa(MatrixXd::Zero(Max+1, Max+1)),tempadag(MatrixXd::Zero(Max+1, Max+1));
                MatrixXd tempeye(MatrixXd::Identity(Max+1, Max+1));

                for(int i=0; i<Max; ++i)
                {
                        tempa(i, i+1)=sqrt(i+1);
                }//cout<<tempa<<endl;
                char n; //cin>>n;
                for(int i=1; i<=Max; ++i)
                {
                        tempadag(i, i-1)=sqrt(i);
                }//cout<<tempadag<<endl;
                //cin>>n;

                sigmaz=kron(tempeye, tempsigmaz);//cout<<sigmaz<<endl;cin>>n;
                sigmamin=kron(tempeye, tempsigmamin);
                sigmaplus=kron(tempeye, tempsigmaplus);

                a=kron(tempa, tempsigmaeye);//cout<<a<<endl;cin>>n;
                adag=kron(tempadag, tempsigmaeye);

                MatrixXd H;

                H=deltar*kron(tempadag*tempa, tempsigmaeye)+deltaq/2*sigmaz
                        +(kron(tempadag, tempsigmamin)+kron(tempa, tempsigmaplus))+f*(adag+a);
                        //cout<<H<<endl;cin>>n;

                

                MatrixXd Iden(kron(tempeye, tempsigmaeye));

                MatrixXcd A;

                A=-i*kron(H, Iden)+i*kron(Iden, H.transpose())+kappa/2*kron(a, adag.transpose())-kappa/2*
                    kron(adag*a, Iden)-kappa/2*kron(Iden, (adag*a).transpose())+gamma*
                    kron(sigmamin, sigmaplus.transpose())-gamma/2*kron(sigmaplus*sigmamin, Iden)-gamma/2*
                    kron(Iden, (sigmaplus*sigmamin).transpose());



                ComplexEigenSolver<MatrixXcd> ces(A);
                //cout<<ces.eigenvalues()<<endl;//cout<<ces.eigenvectors().col(0)<<endl;

                outfile<<deltaq<<"\t"<<endl;

                for(int count=0; count<A.rows(); ++count)
                {
                for(int i=0; i<a.rows(); ++i)
                {
                    for(int j=0; j<a.cols(); ++j)
                    {
                        complex<double> temp=ces.eigenvectors().col(count)(i*a.cols()+j);
                        _rho(i,j)=temp;
                    }
                }

                outfile<<real(ces.eigenvalues()(count))<<"\t"<<abs((_rho*a).trace())<<endl;;

                }
        };

        MatrixXd kron(const MatrixXd& A, const MatrixXd& B)
        {
                MatrixXd kronMatrix(MatrixXd::Zero(A.rows()*B.rows(), A.cols()*B.cols()));
                for(int i=0; i<A.rows(); ++i)
                {
                        for(int j=0; j<A.cols(); ++j)
                        {
                                kronMatrix.block(i*B.rows(), j*B.cols(), B.rows(), B.cols())=A(i, j)*B;
                        }
                }

                return kronMatrix;
        };
        
     

};






#endif
