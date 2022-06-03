//#include <Rcpp.h>
//#include <algorithm>
//#include <cmath>

#include <RcppArmadillo.h>
//[[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
using namespace arma;


////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::uvec bin_drop (int bin_size, int bin_number, arma::rowvec res_std_max){
  
  vec bin_mean(bin_number);
  uvec a;
  for (int i=0;i<bin_number-1;i++){
    a=regspace<arma::uvec>(i*bin_size,(i+1)*bin_size-1);
    bin_mean(i)=mean(res_std_max(a));
  } 
  
  a=regspace<arma::uvec>((bin_number-1)*bin_size,res_std_max.n_cols-1);
  bin_mean(bin_number-1)=mean(res_std_max(a));
  
  uvec drop_bin=find(bin_mean==max(bin_mean));
  //cast uword to int
  if ((int)drop_bin(0)==bin_number-1){
    a=regspace<arma::uvec>(0,(bin_number-1)*bin_size-1);
  }else{
    uvec pos=regspace<arma::uvec>(0,res_std_max.n_cols-1);
    a=join_cols(find((pos< drop_bin[0]*bin_size)),find(pos>((drop_bin[0]+1)*bin_size-1)));
  }
  return (a);                                
}

////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List NMFCPP( const arma::mat f, int loop) { //output svd
    
  //initiate Lagrangian multiplie    
  Function svd = Environment("package:base")["svd"];
    
  arma::mat Lambda(f.n_rows,f.n_cols);
  Lambda.fill(0.);
  arma::mat fitted(f.n_rows,f.n_cols);
 
  int  i,du,dv;    
  List modelsvd = svd(f,du=1,dv=1);
      
  rowvec d=modelsvd[0];  //when assign value from a list, must define the type;
  rowvec u=modelsvd[1];
  rowvec v=modelsvd[2];
     
  fitted=d(0)*u.t()*v;
  arma::vec K(f.n_rows);
     
 //iterate over the NMF-OA process
      
  for (i=0;i<loop;i++){
      Lambda=Lambda-1./sqrt(loop)*(fitted-f);
      Lambda.elem(find(Lambda<0)).zeros();
      modelsvd = svd(f+Lambda,du=1,dv=1);
      rowvec d=modelsvd[0];
      rowvec u=modelsvd[1];
      rowvec v=modelsvd[2];
      fitted=d(0)*u.t()*v;
  }
    return modelsvd;
}


/*///////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List optiNMFCPP_no_base(arma::mat f, arma::vec normFactor, int loop) {
  
  //initialize output list
  
  int num_sample=f.n_rows;
  vec rho=zeros<vec>(num_sample);
  vec K=zeros<vec>(num_sample);
  int convergence=0; // 0 means nothing was done
  List output;
  arma::mat ffilt, res, fitted;
  rowvec E;    
  
  
  //list output initialization
  output("rho")=rho;
  output("convergence")=convergence;
  output("K")=K;
  f.each_col() /= normFactor;
  
  //f is sample by length, each row is a sample
  //filter out bases where the coverage are all low across all samples      
  rowvec colmax = max(f);     //column max
  uvec filter=find(colmax>0.1*colmax.max());
  ffilt=f.cols(filter);    
  arma::uvec nonZero=find(sum(ffilt,1)>0);
  
  //if gene is too short or samples contain 0 read count, return
  //otherwise perform matrix factorization over-approx.
  
  if(filter.size()<50||nonZero.size()<num_sample){
    return output;
  } 
  List modelsvd=NMFCPP(ffilt,loop);
  rowvec d=modelsvd[0];
  rowvec u=modelsvd[1];
  rowvec v=modelsvd[2];
  fitted=d(0)*u.t()*v; 
  uvec q=find(fitted<ffilt);
  if(q.size()>0){         
     fitted(q)=ffilt(q);     
  }
  K=d(0)*abs(u.t());
  rho = 1 - sum(ffilt,1)/(sum(fitted,1) + 1);

  f.each_col() /= K; //each row divide by K    
  E=max(f);  //envelop function
    
  convergence=2 ;  //2 means no-baseline selection;
  output("rho")=rho;
  output("convergence")=convergence;
  output("K")=K;
  output("envelop")=E;    
  return output;
}  
*/       


////////////////////////////////////////////////////////////////////
// optiNMFCPP without binning
// [[Rcpp::export]]
List optiNMFCPP(arma::mat f, arma::vec normFactor, int loop, int baseline) {
        
    //initialize output list
    
    int num_sample=f.n_rows;
    //rho is degradation index
    vec rho=zeros<vec>(num_sample);
    //K is the true scale after degradation normalization
    vec K=zeros<vec>(num_sample);
    //En is envelop function
    rowvec E;
    int convergence=0; // 0 means nothing was done
    List output;
    arma::mat ffilt, res, fitted;
    //list output initialization
    output("rho")=rho;
    output("convergence")=convergence;
    output("K")=K;

    //normalize sequencing depth based on the norm.factor for each sample
    f.each_col() /= normFactor;

    //f is sample by length, each row is a sample
    //filter out bases where the coverage are all low across all samples      
    rowvec colmax = max(f);     //column max
    uvec filter=find(colmax>0.1*colmax.max());
    ffilt=f.cols(filter);    
    int gene_Length=ffilt.n_cols;            
    arma::uvec nonZero=find(sum(ffilt,1)>0);

    //perform non-baseline selection NMF
    // if there are still more than 50 bp left after filtering, proceed; 
    // otherwise return default values
    if(filter.size()<50||(int) nonZero.size()<num_sample){
        return output;
    } 
    
    //perform matrix factorization over approximation on filtered data
    List modelsvd=NMFCPP(ffilt,loop);
    rowvec d=modelsvd[0];
    rowvec u=modelsvd[1];
    rowvec v=modelsvd[2];
    fitted=d(0)*u.t()*v;
    
    //in case some fitted values< raw coverage value, replace fitted by raw value
    
    uvec q=find(fitted<ffilt);
    if(q.size()>0){
        fitted(q)=ffilt(q);
    }
    K=d(0)*abs(u.t());
    fitted=d(0)*u.t()*v; 
    
    res=fitted-ffilt;
    rho = 1 - sum(ffilt,1)/(sum(fitted,1) + 1);

    // if gene_length after filtering is short, or the minimum degradation ratio is >0.2,
    // which means all samples are degraded, do not perform baseline seleciton. 

    if(gene_Length<200||min(rho)>0.2||baseline==0){
        convergence=2 ;  //2 means no-baseline selection;
        output("rho")=rho;
        output("convergence")=convergence;
        output("K")=K;
     }else{
       
        // baseline selection divides genes into 20 bins
        // in each iteration, drop one bin with largest standardized residual

        int bin_size=ceil(gene_Length/20.);
        int bin_number;
        if (gene_Length % bin_size==0){
            bin_number=gene_Length/bin_size;
        }else{
            bin_number=floor(gene_Length/bin_size)+1;
        }


        //a successful baseline is a region that each sample has similar shape
        //or the degradation index score close to 0 for each sample
        
        while(max(rho)>0.1){
            //cout<<"base"<<std::endl;
            arma::mat res_std(res.n_rows,res.n_cols);
            res_std=pow(res/(ffilt+1.0),2);
            rowvec res_std_max=max(res_std,0);
            
            //call function bin_drop to drop bin with largest standardized residual
            uvec a=bin_drop(bin_size, bin_number,res_std_max);
            ffilt=ffilt.cols(a);
            List modelsvd=NMFCPP(ffilt,loop);
            rowvec d=modelsvd[0];
            rowvec u=modelsvd[1];
            rowvec v=modelsvd[2];
            fitted=d(0)*u.t()*v; 
            
            //coersion of over-approximation
            uvec q=find(fitted<ffilt);
            if(q.size()>0){         
               fitted(q)=ffilt(q);     
            }
            K=d(0)*abs(u.t());
            res=fitted-ffilt;
            rho = 1 - sum(ffilt,1)/(sum(fitted,1) + 1);
            bin_number=bin_number-1;
            
            //stop iteration if only there are 6 bin remaining or gene_length is short as 200 
            //or fitted value are close to 0 due to abnormality of algorithm
            if(bin_number==6||ffilt.n_cols<200||min(sum(fitted,1)) < 0.01){
                break;
            }                           
        }
        if(max(rho)<0.2){
          
            //convergence of baseline selection is declared if max of degradation is <0.2
            //use the K at covergence as true K to calculate the enverlop function G
            
            output["convergence"]=1;  //baseline selection successfull
            ffilt=f.cols(filter);  
            rowvec G;
            ffilt.each_col() /= K; //each row divide by K
            G=max(ffilt);  //column max
            fitted=K*G;
            ffilt=f.cols(filter);  
            rho=1-sum(ffilt,1)/(sum(fitted,1)+1);
            output["K"]=K;
            output("rho")=rho; 
            
            //If max of degradation is >0.9, this implies the transcript is extremely degraded
            // and the results are unstable or not reliable, baseline selection results are abandoned
            if(max(rho)>0.9){
               output["convergence"]=3;   //means baseline was found, but degradation too large
            }
          }else{
            output["convergence"]=4;       // baseline selection didn't converge.
          }
         }                

    f.each_col() /= K; //each row divide by K    
    E=max(f);  //envelop function
    output("K")=K;
    output("rho")=rho;
    output("envelop")=E;    
    return output;
}         


////////////////////////////////////////////////////////////////////
    // optiNMFCPP_grid, f is sub-samapled for base line selection
    // grid_size defines the bin size where mean was calculated for baseline
///////////////////////////////////////////////////////////////////    
    // [[Rcpp::export]]
List optiNMFCPP_grid(arma::mat f, arma::vec normFactor, int loop, int grid_size) {
      
       //initialize output list
      int num_sample=f.n_rows;
      vec rho=zeros<vec>(num_sample);
      vec K=zeros<vec>(num_sample);
      rowvec E;
      int convergence=0; // 0 means nothing was done
      List output;
      arma::mat ffilt, res, fitted;
      output("rho")=rho;
      output("convergence")=convergence;
      output("K")=K;
      
      //normalize f by scale normalization factor
      f.each_col() /= normFactor;
      
      //f is sample by length, each row is a sample
      //filter out bases where the coverage are all low across all samples      
      rowvec colmax = max(f);     //column max
      uvec filter=find(colmax>0.1*colmax.max());
      ffilt=f.cols(filter);    
      int gene_Length=ffilt.n_cols;  
      //cout<<"gene_length"<<gene_Length<<std::endl;
      arma::uvec nonZero=find(sum(ffilt,1)>0);

      //perform non-baseline selection NMF
      
      if(filter.size()<50||(int) nonZero.size()<num_sample){
        return output;
      } 
      
      //results to output if without baseline selection
      List modelsvd=NMFCPP(ffilt,loop);
      rowvec d=modelsvd[0];
      rowvec u=modelsvd[1];
      rowvec v=modelsvd[2];
      fitted=d(0)*u.t()*v; 
      uvec q=find(fitted<ffilt);
      if(q.size()>0){         
         fitted(q)=ffilt(q);     
      }
      K=d(0)*abs(u.t());
      res=fitted-ffilt;
      rho = 1 - sum(ffilt,1)/(sum(fitted,1) + 1);
      
      if(gene_Length<200||min(rho)>0.2){
        convergence=2 ;  //2 means no-baseline selection;
        output("rho")=rho;
        output("convergence")=convergence;
        output("K")=K;
      }else{
        
        //baseline selection after on grids
        int grid_number;
        if(gene_Length/grid_size<20){
          grid_size = floor(gene_Length/20);
        }
        //binning fcolumns to improve efficiency
        if (gene_Length % grid_size==0){
          grid_number=gene_Length/grid_size;
        }else{
          grid_number=floor(gene_Length/grid_size)+1;
        }
          
        arma::mat ffilt_bin(ffilt.n_rows,grid_number);
        //cout<<"grid_number"<<grid_number<<std::endl;
        
        uvec a;
        for (int i=0;i<grid_number-1;i++){
          a=regspace<arma::uvec>(i*grid_size,(i+1)*grid_size-1);
          ffilt_bin.col(i)=mean(ffilt.cols(a),1);
         } 
        
        a=regspace<arma::uvec>((grid_number-1)*grid_size,ffilt.n_cols-1);
        ffilt_bin.col(grid_number-1)=mean(ffilt.cols(a),1);
        //cout<<"here"<<std::endl;   
        
        List modelsvd=NMFCPP(ffilt_bin,loop);
        //cout<<"here"<<std::endl;
        
        rowvec d=modelsvd[0];
        rowvec u=modelsvd[1];
        rowvec v=modelsvd[2];
        fitted=d(0)*u.t()*v; 
        K=d(0)*abs(u.t());
        uvec q=find(fitted<ffilt_bin);
        if(q.size()>0){         
           fitted(q)=ffilt_bin(q);     
        }
        //cout<<"here"<<std::endl;
        
        res=fitted-ffilt_bin;
        rho = 1 - sum(ffilt_bin,1)/(sum(fitted,1) + 1);
        
        //enter baseline selection, splitting ffilt_bin into ~20 bins
        
        int bin_size=ceil(ffilt_bin.n_cols/20.);
        int bin_number;
        if (ffilt_bin.n_cols % bin_size==0){
          bin_number=ffilt_bin.n_cols/bin_size;
        }else{
          bin_number=floor(ffilt_bin.n_cols/bin_size)+1;
        }
        
        while(max(rho)>0.1){
          //cout<<"base"<<std::endl;
          arma::mat res_std(res.n_rows,res.n_cols);
          res_std=pow(res/(ffilt_bin+1.0),2);
          rowvec res_std_max=max(res_std,0);
           
          //call function bin_drop to drop bin with largest standardized residual
          uvec a=bin_drop(bin_size, bin_number,res_std_max);
          //cout<<"here"<<std::endl;
          
          ffilt_bin=ffilt_bin.cols(a);
          List modelsvd=NMFCPP(ffilt_bin,loop);

          rowvec d=modelsvd[0];
          rowvec u=modelsvd[1];
          rowvec v=modelsvd[2];
          fitted=d(0)*u.t()*v; 

          //coersion of over-approximation
          uvec q=find(fitted<ffilt_bin);
          if(q.size()>0){         
            fitted(q)=ffilt_bin(q);     
          }
          K=d(0)*abs(u.t());
          res=fitted-ffilt_bin;
          rho = 1 - sum(ffilt_bin,1)/(sum(fitted,1) + 1);
          bin_number=bin_number-1;

          if(bin_number<0.3*gene_Length/grid_size/bin_size||
             (int) ffilt_bin.n_cols<200/grid_size/bin_size||
             min(sum(fitted,1)) < 0.01){
            break;
          }                           
        }
        if(max(rho)<0.2){
          output["convergence"]=1;  //baseline selection successfull
          
          ffilt_bin=f.cols(filter);  
          rowvec G;
          ffilt.each_col() /= K; //each row divide by K
          G=max(ffilt);  //column max
          fitted=K*G;
          ffilt=f.cols(filter);  
          rho=1-sum(ffilt,1)/(sum(fitted,1)+1);
          
          output["K"]=K;
          output("rho")=rho; 
          
          if(max(rho)>0.9){
            output["convergence"]=3;   //means baseline was found, but degradation too large
          }
        }else{
          output["convergence"]=4;       // baseline selection didn't converge.
        }
      }                
      f.each_col() /= K; //each row divide by K    
      E=max(f);  //envelop function
      output("K")=K;
      output("rho")=rho;
      output("envelop")=E;     
      return output;
}         


