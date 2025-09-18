module constants
implicit none
integer dim,neln,nparam,intn,nnode,nel,nframe,di,group,npop,lam_length
integer sqnn, break
real(8) tweight,alpha
real(8) lowb(14),upb(14),s_sta,c_sta_k1(5),c_sta_k2(5)
real(8) xrate,mrate,xcon,mcon,keep_prob,Tmax,Tmin
integer nbit,ngen,nt,mnum,nstac
integer mut_mode,cross_mode,uni_penalty
integer control_boltz,elitism
integer phys_c_e,phys_c_l,bctype,run_mode
integer ncomp,nten
integer rflag,nnode2,nel2
parameter(dim=2,neln=4,nparam=14,intn=1,nnode=416,nel=375)
parameter(nnode2=390,nel2=350,ncomp=17,nten=13,break=18)
parameter(nframe=30,di=8,group=4,npop=100,sqnn=52)
parameter(alpha=1.400d0,nbit=4,ngen=80,nt=300)
parameter(xrate=0.840d0,mrate=0.030d0,xcon=0.900d0,mcon=0.200d0,mnum=1)
parameter(nstac = 5, s_sta = 0.510d0,tweight = 0.700d0)
parameter(c_sta_k1=(/-0.150d0,-0.250d0,-0.350d0,-0.450d0,-0.550d0/))
parameter(c_sta_k2=(/0.150d0,0.250d0,0.350d0,0.450d0,0.550d0/))
parameter(mut_mode = 1 ,cross_mode = 1,uni_penalty = 0)
parameter(control_boltz = 1 ,phys_c_e = 1, phys_c_l = 3)
parameter(Tmax = 60.000d0,Tmin = 0.500d0)
parameter(elitism = 1)
parameter(bctype = 2)
parameter(run_mode = 2)
parameter(rflag = 1)
parameter(lowb=(/82000.000d0,153800.000d0,0.100d0,0.050d0,-0.000d0,0.050d0,2.000d0,0.010d0,0.050d0, &
 &    2.000d0,2.000d0,5.000d0,0.001d0,1.000d0/))
parameter(upb=(/122000.000d0,233800.000d0,0.300d0,6.000d0,-0.500d0,0.400d0,10.000d0,1.000d0,6.000d0, &
 &    8.000d0,8.000d0,25.000d0,0.500d0,6.000d0/))
end module constants