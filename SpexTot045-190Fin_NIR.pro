function Fin,X,P

  wnorm = 8
  ocdir='/Users/mmarsset/IDLWorkspace85/Default/Shkuratov/ocs/'
  
  
  openr,1,ocdir+'oliv70.dat'
  data_ol=fltarr(3,59)
  readf,1,data_ol
  close,1
  
  data_ol=data_ol(*,15:58)
  
  openr,2,ocdir+'opx70.dat'
  data_opx=fltarr(3,59)
  readf,2,data_opx
  close,2
  
  data_opx=data_opx(*,15:58)
  
  openr,3,ocdir+'cpx70.dat'
  data_cpx=fltarr(3,59)
  readf,3,data_cpx
  close,3
  
  data_cpx=data_cpx(*,15:58)
  
  openr,5,ocdir+'iron.dat'
  data5=fltarr(3,59)
  readf,5,data5
  close,5
  
  data5=data5(*,15:58)
  
  openr,6,ocdir+'chromite.dat'
  data6=fltarr(3,59)
  readf,6,data6
  close,6
  
  data6=data6(*,15:58)
  
  
  n_ol=fltarr(44)
  n_ol=data_ol(1,*)
  n_opx=fltarr(44)
  n_opx=data_opx(1,*)
  n_cpx=fltarr(44)
  n_cpx=data_cpx(1,*)
  n_ir=fltarr(44)
  n_ir=data5(1,*)
  n_irOx=fltarr(44)
  n_irOx=data6(1,*)
  
  ro_ir=fltarr(44)
  Re_ir=fltarr(44)
  Ri_ir=fltarr(44)
  Rb_ir=fltarr(44)
  Rf_ir=fltarr(44)
  Te_ir=fltarr(44)
  Ti_ir=fltarr(44)
  
  ro_irOx=fltarr(44)
  Re_irOx=fltarr(44)
  Ri_irOx=fltarr(44)
  Rb_irOx=fltarr(44)
  Rf_irOx=fltarr(44)
  Te_irOx=fltarr(44)
  Ti_irOx=fltarr(44)
  
  ro_ol=(n_ol-1)^2/(n_ol+1)^2
  Re_ol=ro_ol+0.05
  Ri_ol=1.04-1/(n_ol^2)
  Rb_ol=(0.28*n_ol-0.2)*Re_ol
  Rf_ol=Re_ol-Rb_ol
  Te_ol=1-Re_ol
  Ti_ol=Te_ol/(n_ol^2)
  
  ro_opx=(n_opx-1)^2/(n_opx+1)^2
  Re_opx=ro_opx+0.05
  Ri_opx=1.04-1/(n_opx^2)
  Rb_opx=(0.28*n_opx-0.2)*Re_opx
  Rf_opx=Re_opx-Rb_opx
  Te_opx=1-Re_opx
  Ti_opx=Te_opx/(n_opx^2)
  
  ro_cpx=(n_cpx-1)^2/(n_cpx+1)^2
  Re_cpx=ro_cpx+0.05
  Ri_cpx=1.04-1/(n_cpx^2)
  Rb_cpx=(0.28*n_cpx-0.2)*Re_cpx
  Rf_cpx=Re_cpx-Rb_cpx
  Te_cpx=1-Re_cpx
  Ti_cpx=Te_cpx/(n_cpx^2)
  
  t_ol=fltarr(1,44)
  t_opx=fltarr(1,44)
  t_cpx=fltarr(1,44)
  t_ir=fltarr(1,44)
  t_irOx=fltarr(1,44)
  r_ol_b=fltarr(1,44)
  r_ol_f=fltarr(1,44)
  r_opx_b=fltarr(1,44)
  r_opx_f=fltarr(1,44)
  r_cpx_b=fltarr(1,44)
  r_cpx_f=fltarr(1,44)
  r_ir_b=fltarr(1,44)
  r_ir_f=fltarr(1,44)
  r_irOx_b=fltarr(1,44)
  r_irOx_f=fltarr(1,44)
  pb=fltarr(1,44)
  pf=fltarr(1,44)
  A=fltarr(44)
  F=fltarr(1,44)
  
  
  Eps_ol=fltarr(2,44)
  Eps_opx=fltarr(2,44)
  Eps_cpx=fltarr(2,44)
  Eps_ir=fltarr(2,44)
  Eps_ol=complex(data_ol(1,*),data_ol(2,*),/Double)
  Eps_opx=complex(data_opx(1,*),data_opx(2,*),/Double)
  Eps_cpx=complex(data_cpx(1,*),data_cpx(2,*),/Double)
  Eps_ir=complex(data5(1,*),data5(2,*),/Double)
  
  Div_ol=fltarr(2,44)
  Div_opx=fltarr(2,44)
  Div_cpx=fltarr(2,44)
  k_ol_weath=fltarr(1,44)
  k_opx_weath=fltarr(1,44)
  k_cpx_weath=fltarr(1,44)
  Div_ol=[(Eps_ir(*)^2/Eps_ol(*)^2)-1]/[(Eps_ir(*)^2/Eps_ol(*)^2)+2]
  Div_opx=[(Eps_ir(*)^2/Eps_opx(*)^2)-1]/[(Eps_ir(*)^2/Eps_opx(*)^2)+2]
  Div_cpx=[(Eps_ir(*)^2/Eps_cpx(*)^2)-1]/[(Eps_ir(*)^2/Eps_cpx(*)^2)+2]
  
  q=0.5
  
  ;k_ol_weath=(3*P(6)*n_ol/20)*imaginary(Div_ol(*))
  ;k_opx_weath=(3*P(6)*n_opx/20)*imaginary(Div_opx(*))
  ;k_cpx_weath=(3*P(6)*n_cpx/20)*imaginary(Div_cpx(*))
  
  ro_ir(*)=(n_ir(*)-1)^2/(n_ir(*)+1)^2
  Re_ir(*)=ro_ir(*)+0.05
  Ri_ir(*)=1.04-1/(n_ir(*)^2)
  ;Rb_ir(*)=(0.28*n_ir(*)-0.2)*Re_ir(*)/2.25+0.1
  Rb_ir(*)=(0.28*n_ir(*)-0.2)*(Re_ir(*))/2.3+0.22
  Rf_ir(*)=Re_ir(*)-Rb_ir(*)
  Te_ir(*)=1-Re_ir(*)
  Ti_ir(*)=Te_ir(*)/(n_ir(*)^2)
  
  
  ro_irOx(*)=(n_irOx(*)-1)^2/(n_irOx(*)+1)^2
  Re_irOx(*)=ro_irOx(*)+0.05
  Ri_irOx(*)=1.04-1/(n_irOx(*)^2)
  Rb_irOx(*)=(0.28*n_irOx(*)-0.2)*Re_irOx(*)
  ;Rb_irOx(*)=(0.28*n_irOx(*)-0.2)*Re_irOx(*)/1.95+0.17 ; Nickel
  Rf_irOx(*)=Re_irOx(*)-Rb_irOx(*)
  Te_irOx(*)=1-Re_irOx(*)
  Ti_irOx(*)=Te_irOx(*)/(n_irOx(*)^2)
  
  ;t_ol(*)=(4*!Pi*(data_ol(2,*)+k_ol_weath(*))*P(5))/data_ol(0,*)
  t_ol(*)=(4*!Pi*data_ol(2,*)*P(5))/data_ol(0,*)
  ;t_opx(*)=(4*!Pi*(data_opx(2,*)+k_opx_weath(*))*P(5))/data_opx(0,*)
  t_opx(*)=(4*!Pi*data_opx(2,*)*P(5))/data_opx(0,*)
  t_cpx(*)=(4*!Pi*data_cpx(2,*)*P(5))/data_cpx(0,*)
  t_ir(*)=(4*!Pi*data5(2,*)*P(5))/data5(0,*)
  t_irOx(*)=(4*!Pi*data6(2,*)*P(5))/data6(0,*)
  
  r_ol_b(*)=Rb_ol+0.5*Te_ol*Ti_ol*Ri_ol*exp(-2*t_ol(*))/(1-Ri_ol*exp(-t_ol(*)))
  r_ol_f(*)=Rf_ol+Te_ol*Ti_ol*exp(-t_ol(*))+0.5*Te_ol*Ti_ol*Ri_ol*exp(-2*t_ol(*))/(1-Ri_ol*exp(-t_ol(*)))
  r_opx_b(*)=Rb_opx+0.5*Te_opx*Ti_opx*Ri_opx*exp(-2*t_opx(*))/(1-Ri_opx*exp(-t_opx(*)))
  r_opx_f(*)=Rf_opx+Te_opx*Ti_opx*exp(-t_opx(*))+0.5*Te_opx*Ti_opx*Ri_opx*exp(-2*t_opx(*))/(1-Ri_opx*exp(-t_opx(*)))
  r_cpx_b(*)=Rb_cpx+0.5*Te_cpx*Ti_cpx*Ri_cpx*exp(-2*t_cpx(*))/(1-Ri_cpx*exp(-t_cpx(*)))
  r_cpx_f(*)=Rf_cpx+Te_cpx*Ti_cpx*exp(-t_cpx(*))+0.5*Te_cpx*Ti_cpx*Ri_cpx*exp(-2*t_cpx(*))/(1-Ri_cpx*exp(-t_cpx(*)))
  r_ir_b(*)=Rb_ir(*)+0.5*Te_ir(*)*Ti_ir(*)*Ri_ir(*)*exp(-2*t_ir(*))/(1-Ri_ir(*)*exp(-t_ir(*)))
  r_ir_f(*)=Rf_ir(*)+Te_ir(*)*Ti_ir(*)*exp(-t_ir(*))+0.5*Te_ir(*)*Ti_ir(*)*Ri_ir(*)*exp(-2*t_ir(*))/(1-Ri_ir(*)*exp(-t_ir(*)))
  r_irOx_b(*)=Rb_irOx(*)+0.5*Te_irOx(*)*Ti_irOx(*)*Ri_irOx(*)*exp(-2*t_irOx(*))/(1-Ri_irOx(*)*exp(-t_irOx(*)))
  r_irOx_f(*)=Rf_irOx(*)+Te_irOx(*)*Ti_irOx(*)*exp(-t_irOx(*))+0.5*Te_irOx(*)*Ti_irOx(*)*Ri_irOx(*)*exp(-2*t_irOx(*))/(1-Ri_irOx(*)*exp(-t_irOx(*)))
  
  
  pb(*)=q*(P(0)*r_ol_b(*)+P(1)*r_opx_b(*)+P(2)*r_cpx_b(*)+P(3)*r_ir_b(*)+P(4)*r_irOx_b(*))/(P(0)+P(1)+P(2)+P(3)+P(4))
  pf(*)=q*(P(0)*r_ol_f(*)+P(1)*r_opx_f(*)+P(2)*r_cpx_f(*)+P(3)*r_ir_f(*)+P(4)*r_irOx_f(*))/(P(0)+P(1)+P(2)+P(3)+P(4))+1-q
  
  F(*)=exp(-P(6)/data_ol(0,*))*[(1+pb(*)^2-pf(*)^2)/(2*pb(*))-[((1+pb(*)^2-pf(*)^2)/(2*pb(*)))^2-1]^0.5]
  
  
  ;A(*)=transpose(F(*))/F(wnorm)
  A(*)=transpose(F(*))/F(P(7))
  ;A(*)=transpose(F(*))
  return,A

end


;just to get x values for interpolation of asteroid data
ocdir='/Users/mmarsset/IDLWorkspace85/Default/Shkuratov/ocs/'
openr,1,ocdir+'oliv70.dat'
data_ol=fltarr(3,59)
readf,1,data_ol
close,1
data_ol=data_ol(*,15:58)
lenarr=n_elements(data_ol(0,*))


;INPUT SPECTRUM
;**************INPUTS********************
dir = '/Users/mmarsset/data/IRTF/MBA_families_spectra/502_Eunomia/'
fileout=dir+'out_nir_binzel.dat'
ncol=2
;****************************************


readcol,dir+'list_nir.dat',file, format='(A)'

outinfo=fltarr(7,n_elements(file))

for z=0,n_elements(file)-1 do begin

  wnorm = 3

  ;print,'************ '+file(z)
  if ncol eq 2 then readcol,dir+file(z),xarr,yarr,/silent
  if ncol eq 3 then readcol,dir+file(z),xarr,yarr,earr,/silent
  if ncol eq 4 then readcol,dir+file(z),xarr,yarr,earr,foo,/silent
  good=where(yarr gt 0.) ; ignore negatives
  xarr=xarr(good)
  yarr=yarr(good)
  data11=fltarr(2,n_elements(xarr))
  data11(0,*)=xarr
  data11(1,*)=yarr
  
  data10=fltarr(2,lenarr)
  data10(0,*) = data_ol(0,*)
  data10(1,*) = interpol(data11(1,*), data11(0,*),data_ol(0,*))
  
  d=fltarr(lenarr)
  d=transpose(data10(1,*)/data10(1,wnorm))
  ;d=transpose(data10(1,*)*0.25)
  result=d
  
  X=fltarr(lenarr)
  X=transpose(data10(0,*))
  Y=fltarr(lenarr)
  Y=result
  P=[1.,1,0.0,0.0,0.1,10,0.1,wnorm]
  weights=1./Y
  
  parinfo=replicate({value:0,fixed:0,limited:[1,0],limits:[0.D,1.D]},8)
  
  parinfo(0).fixed=0
  parinfo(1).fixed=0
  parinfo(2).fixed=1
  parinfo(3).fixed=1
  parinfo(4).fixed=1
  parinfo(5).fixed=0
  parinfo(6).fixed=0
  parinfo(7).fixed=1
  
  ;parinfo(0).limited(0)=1
  parinfo(0).limits(0)=0.D
  parinfo(0).limits(1)=1.D
  
  ;parinfo(1).limited(0)=1
  parinfo(1).limits(0)=0.D
  parinfo(1).limits(1)=1.D
  
  ;parinfo(2).limited(0)=1
  parinfo(2).limits(0)=0.D
  parinfo(2).limits(1)=1.D
  
  ;parinfo(3).limited(0)=1
  parinfo(3).limits(0)=0.D
  parinfo(3).limits(1)=1.D
  
  ;parinfo(4).limited(0)=1
  parinfo(4).limits(0)=0.D
  parinfo(4).limits(1)=1.D
  
  ;parinfo(5).limited(0)=1
  parinfo(5).limits(0)=1.D
  parinfo(5).limits(1)=300.D
  
  ;parinfo(6).limited(0)=1
  parinfo(6).limits(0)=-1.D
  parinfo(6).limits(1)=1.D
  
  parinfo(7).limits(0)=0
  parinfo(7).limits(1)=100  
  
  
  yfit=mpfitfun('Fin',X,Y,weights,P,PARINFO=parinfo)
  
  outinfo(0,z)=70
  outinfo(1,z)=yfit(0)/(yfit(0)+yfit(1))
  outinfo(2,z)=yfit(0)
  outinfo(3,z)=yfit(1)
  outinfo(4,z)=0.1
  outinfo(5,z)=yfit(5)
  outinfo(6,z)=yfit(6)
  
  A=fin(data_ol(0,*),yfit)
  norfac = mean(A(*))/mean(data10(1,*))
  p1 = plot(data10(0,*),data10(1,*),xrange=[0.4,1.95],/xstyle,yrange=[min(data10(1,*))-0.2,max(data10(1,*))+0.2],/ystyle)
  p2 = plot(data10(0,*),A(*)/norfac,/xstyle,/ystyle,color='red',/overplot,xtitle='Wavelength [micron]',ytitle='Reflectance',title=file[z])
  ;p2.save, dir+'/'+file[z].replace('.txt','.png')

endfor

fmt='A40,7F10.4'
cmt='         file             Mg    ol/ol+opx     olv      opx       chr       GS        CS '
forprint,file,outinfo(0,*),outinfo(1,*),outinfo(2,*),outinfo(3,*),outinfo(4,*),outinfo(5,*),outinfo(6,*),format=fmt,textout=fileout,comment=cmt



end
