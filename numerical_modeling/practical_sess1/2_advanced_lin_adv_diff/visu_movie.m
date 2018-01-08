function M=visu_movie(phi,numfig)
  nt=size(phi,1);
  for it=1:nt
  figure(numfig)
  %clf(numfig)
  plot(phi(it,:),'b');
  set(gca,'xlim',[1 size(phi,2)]);
  set(gca,'ylim',[-0.05 1.5]);
  M(it)=getframe;
  end   
 %movie(M,1);
