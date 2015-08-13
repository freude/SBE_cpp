/*
 * Heterostr.cpp
 *
 *  Created on: Mar 12, 2012
 *      Author: mk
 */

#include "Heterostr.h"
#include "TernaryAlloy.h"
#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <iostream>
#include <cmath>
#include "my_constants.h"
#include "my_functions.h"


namespace std {

Heterostr::Heterostr() {
	int flag=0;
	Temper=300;
	heterostr="AlAs/GaAs/AlAs";
	NumLayers=3;
	rat=0.3;
	widtharray = new double[NumLayers];
	for (size_t j=0;j<NumLayers;j++) widtharray[j]=5e-9;
	flag=Het2Materials(heterostr);
}

Heterostr::Heterostr(string het, Vec_DP& width) {
	int flag=0;
	Temper=300;
	heterostr=het;
	NumLayers=width.size();
	rat=0.3;
	widtharray = new double[NumLayers];
	for (size_t j=0;j<NumLayers;j++) widtharray[j]=width(j);
	flag=Het2Materials(heterostr);
}

Heterostr::Heterostr(string het, Vec_DP& width, double T) {
	int flag=0;
	Temper=T;
	heterostr=het;
	NumLayers=width.size();
	rat=0.3;
	widtharray = new double[NumLayers];
	for (size_t j=0;j<NumLayers;j++) widtharray[j]=width(j);
	flag=Het2Materials(heterostr);
}

Heterostr::~Heterostr() {
	delete [] matarray;
	delete [] widtharray;

}

int Heterostr::Het2Materials(const string str){

	int flag=1;
	size_t pos=0, pos1=0, j1=0;;

	matarray = new TernaryAlloy[NumLayers];

	pos=str.find("/");

	if ((pos==string::npos)||(NumLayers<3)){
		flag=0;
	}
	else{
		while (pos!=string::npos){
			j1=j1+1;
			matarray[j1-1].IntParams(str.substr(pos1,pos-pos1),Temper);
			pos1=pos+1;
			pos=str.find("/",pos1);
		}
		matarray[j1].IntParams(str.substr(pos1),Temper);
		if (j1==NumLayers){
			flag=1;
		}
	}
	return flag;
}

void Heterostr::PrintStr(){
	cout << "---------------------------------"<<endl;
	cout << "The heterostructure is consist of"<<endl;
	cout << "---------------------------------"<<endl;
	cout << "   Material      Width (nm)   "<<endl;
	cout << "---------------------------------"<<endl;
	for (size_t j=0;j<NumLayers;j++){
		cout<<"|    "<<matarray[j].GetMaterial()<<"   |   "<<widtharray[j]<<"    |"<<endl;
	}
	cout << "---------------------------------"<<endl;
	return;
}

Vec_DP Heterostr::getParam(Vec_DP& x,string par){
	double q=1.6e-19;
	double m0=9.11e-31;
	size_t l_x=x.size();
	Vec_DP temp(NumLayers);
	Vec_DP param(l_x);

	if (strcmp("a",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("a");
		}
		param=1e-10*pos_dep(x,temp);
	}
	else if (strcmp("Eg",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("Eg");
		}
		param=q*pos_dep(x,temp);
	}
	else if (strcmp("me",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("me");
		}
		param=m0*pos_dep(x,temp);
	}
	else if (strcmp("gamma1",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("gamma1");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("gamma2",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("gamma2");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("gamma3",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("gamma3");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("Ep",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("Ep");
		}
		param=q*pos_dep(x,temp);
	}
	else if (strcmp("delta",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("delta");
		}
		param=q*pos_dep(x,temp);
	}
	else if (strcmp("ac",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("ac");
		}
		param=q*pos_dep(x,temp);
	}
	else if (strcmp("av",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("av");
		}
		param=q*pos_dep(x,temp);
	}
	else if (strcmp("b",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("b");
		}
		param=q*pos_dep(x,temp);
	}
	else if (strcmp("d",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("d");
		}
		param=q*pos_dep(x,temp);
	}
	else if (strcmp("C11",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("C11");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("C12",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("C12");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("C44",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("C44");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("alpha",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("alpha");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("betha",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("betha");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("VBO",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("VBO");
		}
		param=pos_dep(x,temp);
	}
	else if (strcmp("Ec",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("Eg");
		}
		param=q*pos_dep(x,temp);
		//p=rat*(Eg-min(Eg))*q+q*min(Eg);
	}
	else if (strcmp("Ev",par.c_str())==0){
		for (size_t j=0;j<NumLayers;j++){
			temp(j)=matarray[j].GetParam("Eg");
		}
		param=q*pos_dep(x,temp);
		//p=-(1-rat)*(Eg-min(Eg))*q;
	}

	else cout<<"Parameter "<<par<<" is absent in the database";

	return param;
}

Vec_DP Heterostr::pos_dep(Vec_DP& x, Vec_DP& par){

    double spread=1e12;
	Vec_DP Eva(x.size(),par[0]);
    double xb=widtharray[0];

    for (size_t j=1;j<(NumLayers);j++){
    	Eva=Eva+(par[j]-par[j-1])/2*tanh(spread*(x-xb))+(par[j]-par[j-1])/2;
        xb=xb+widtharray[j];
    }
    return Eva;
}

void Heterostr::computeBS(BS& bs, string typ, Compar par){

	if (strcmp("4x4pw",typ.c_str())==0){
		void qwzb4x4ps(BS& bs, Compar par);
	}
	else if (strcmp("4x4pw",typ.c_str())==0){
//		void qwzb6x6ps(BS& bs, string typ);
	}
	else cout<<"Parameter "<<typ<<" is absent in the database";

}

void Heterostr::qwzb4x4ps(BS& bs, Compar par){

	int flag;

	Vec_DP Ev(bs.x.size());
	Vec_DP gamma1(bs.x.size());
	Vec_DP gamma2(bs.x.size());
	Vec_DP gamma3(bs.x.size());

	Ev=getParam(bs.x,"Vh");
	gamma1=getParam(bs.x,"gamma1");
	gamma2=getParam(bs.x,"gamma2");
	gamma3=getParam(bs.x,"gamma3");

	//-------------------------------------------------------------------
	//-----------------------Forier transformations----------------------
	//-------------------------------------------------------------------

	double NN=50;

	double L=bs.x[bs.x.size()-1]-bs.x[0];
	double k=2*pi/L;

	Vec_DP n(-NN,NN,1);
	Vec_DP m(-NN,NN,1);

	Mat_CPLX_DP gamma1f(n.size());
	Mat_CPLX_DP gamma2f(n.size());
	Mat_CPLX_DP gamma3f(n.size());
	Mat_CPLX_DP Evf(n.size());

	double omega;

	for (size_t j1=1;j1<n.size();j1++){
	    for (size_t j2=1;j2<m.size();j2++){

	        omega=k*(m(j2)-n(j1));

	        gamma1f(j1,j2)=trapz(bs.x,gamma1*exp(ij*omega*bs.x))/L;
	        gamma2f(j1,j2)=trapz(bs.x,gamma2*exp(ij*omega*bs.x))/L;
	        gamma3f(j1,j2)=trapz(bs.x,gamma3*exp(ij*omega*bs.x))/L;

	        Evf(j1,j2)=trapz(bs.x,Ev*exp(ij*omega*bs.x))/L;

	    }
	}

	Mat_CPLX_DP M(4*n.size());
	Mat_CPLX_DP FFF[bs.kx.vector->size][bs.ky.vector->size];
	Vec_DP E[bs.kx.vector->size][bs.ky.vector->size];

	gsl_complex P;
	gsl_complex Q;
	gsl_complex S;
	gsl_complex Sc;
	gsl_complex R;
	gsl_complex Rc;

	MatrixCP H2(4);
	MatrixCP H1(4);
	MatrixCP H0(4);
	MatrixCP T(4);

	gsl_eigen_herm_workspace* wksp;

	for (int jj=0;jj<bs.kx.vector->size;jj++){
	    for (int jjj=0;jjj<bs.ky.vector->size;jjj++){


	        for (int j1=0;j1<n.vector->size;j1++){
	        	for (int j2=0;j2<m.vector->size;j2++){

	                P=(pow(h,2))/(2*m0)*gamma1f(j1,j2);
	                Q=-2*(pow(h,2))/(2*m0)*gamma2f(j1,j2);

	                H2(0,0)=-P+Q;
	                H2(1,1)=-P-Q;
	                H2(2,2)=-P-Q;
	                H2(3,3)=-P+Q;

	                S=-ij*sqrt(3).*gamma3f(j1,j2).*(h^2)/m0*(bs.kx(jjj)-ij*bs.ky(jj));
	                Sc=ij*sqrt(3).*gamma3f(j1,j2).*(h^2)/m0*(bs.kx(jjj)+ij*bs.ky(jj));

	                H1(0,1)=-Sc;
	                H1(1,0)=-S;
	                H1(2,3)=Sc;
	                H1(3,2)=S;

	                P=-Evf(j1,j2)+gamma1f(j1,j2)*(h^2)/(2*m0)*(bs.kx(jjj)^2+bs.ky(jj)^2);
	                Q=gamma2f(j1,j2)*(h^2)/(2*m0)*(bs.kx(jjj)^2+bs.ky(jj)^2);
	                R=-sqrt(3)*(h^2)/(2*m0).*(gamma2f(j1,j2)*(bs.kx(jjj)^2-bs.ky(jj)^2)-2*ij*gamma3f(j1,j2)*bs.kx(jjj)*bs.ky(jj));
	                Rc=-sqrt(3)*(h^2)/(2*m0).*(gamma2f(j1,j2)*(bs.kx(jjj)^2-bs.ky(jj)^2)+2*ij*gamma3f(j1,j2)*bs.kx(jjj)*bs.ky(jj));

	                H0(0,0)=-P+Q;
	                H0(0,2)=R;
	                H0(1,1)=-P-Q;
	                H0(1,3)=R;
	                H0(2,0)=Rc;
	                H0(2,2)=-P-Q;
	                H0(3,1)=Rc;
	                H0(3,3)=-P+Q;

	                T=0.5*k*k*(n(j1)*H2*m(j2)+m(j2)*H2*n(j1))+0.5*k*(H1.*m(j2)+n(j1)*H1)+H0;

	                M.assign2submat(4*j1-3, 4*j2-3,T);

	        	}
	        }

	        wksp=gsl_eigen_herm_alloc(4*n.vector->size);
	        flag=gsl_eigen_herm(M.matrix,E[jj][jjj].vector,wksp);
	        gsl_eigen_herm_free(wksp);

	        %     SS=sparse(M);
	        %     [FF,EEE]=eigs(SS,num,'lm');

	        [E(jj,jjj,:), IX]=sort((diag(EEE./q)),1);

	        if cWF_flag>0
	            FFF(jj,jjj,:,:)=FF(:,IX);
	        end;
		}
	}

	//------------------Post-proccessing-------------------------

	[k0,kx0]=min(abs(kx));
	[k0,ky0]=min(abs(ky));

	if only_conf_flag==1
	    index=1;
	    while E(kx0,ky0,end-index+1)>min(Ev/q)
	        index=index+1;
	    end;
	    index=index-1;
	    Eout=squeeze(E(:,:,end-index+1:end));

	    if cWF_flag>0
	        WF=zeros(length(ky),length(kx),index,length(x));
	        for jj=1:length(ky)
	            for jjj=1:length(kx)
	                for j=1:index
	                    for j1=1:length(m)
	                        F1(1,1,1,:)=FFF(jj,jjj,4*j1-3,end-j+1)*exp(ij*k*m(j1)*x);
	                        F2(1,1,1,:)=FFF(jj,jjj,4*j1-2,end-j+1)*exp(ij*k*m(j1)*x);
	                        F3(1,1,1,:)=FFF(jj,jjj,4*j1-1,end-j+1)*exp(ij*k*m(j1)*x);
	                        F4(1,1,1,:)=FFF(jj,jjj,4*j1,end-j+1)*exp(ij*k*m(j1)*x);
	                        WF(jj,jjj,j,:)=WF(jj,jjj,j,:)+F1+F2+F3+F4;
						}
					}
				}
			}
	        varargout(2)=x;
	        varargout(1)=squeeze(WF);
	    }

	else
	    index=SB_num;
	    Eout=squeeze(E(:,:,end-index+1:end));

	    if cWF_flag>0
	        WF=zeros(length(ky),length(kx),index,length(x));
	        for jj=1:length(ky)
	            for jjj=1:length(kx)
	                for j=1:index
	                    for j1=1:length(m)
	                        F1(1,1,1,:)=FFF(jj,jjj,4*j1-3,end-j+1)*exp(ij*k*m(j1)*x);
	                        F2(1,1,1,:)=FFF(jj,jjj,4*j1-2,end-j+1)*exp(ij*k*m(j1)*x);
	                        F3(1,1,1,:)=FFF(jj,jjj,4*j1-1,end-j+1)*exp(ij*k*m(j1)*x);
	                        F4(1,1,1,:)=FFF(jj,jjj,4*j1,end-j+1)*exp(ij*k*m(j1)*x);
	                        WF(jj,jjj,j,:)=WF(jj,jjj,j,:)+F1+F2+F3+F4;
	                    end;
	                end;
	            end;
	        end;
	        varargout(2)=x;
	        varargout(1)=squeeze(WF);
	    end;

	end;

	if save_flag==1
	    fpath=strcat(fpath,'/E.dat');
	    dlmwrite(fpath,Eout,' ');
	    if cWF_flag>0
	        fpath=strcat(fpath,'/WF.dat');
	        dlmwrite(fpath,squeeze(WF(kx0,ky0,:,:)),' ');
	    end;
	end;
}
} /* namespace std */
