function st = machplotter(l_r, x_ob, m_ob, m_fb,q)


            %function to plot pressure after each shock. 
            
            v = cumsum(l_r);
            
                for i = 1:length(v)
                    
                    if(any(v(i) > x_ob(1)))
                        
                        v(i) = x_ob(1);
                    end
                    
                    
                end
                k = [0,v];
                    
                    
                for j = 1 : length(k)-1   
                    plot([k(j:j+1)], m_fb(j)*ones(2,1), 'r-o');
                    hold on;
                end
                
                plot(x_ob(1:q+1), m_ob(1:q+1),'-o');   
                
                hold off;
                
                xlabel('Length from nose to the isolator inlet');
                ylabel('Mach Number Variation');
                legend('ramp','cowl');
              st=0;  
end
                                      