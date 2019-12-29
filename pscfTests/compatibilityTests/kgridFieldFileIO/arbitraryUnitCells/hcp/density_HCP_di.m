function [rho_A, rho_B, rho_C, wave] = density_HCP_di



file_id = fopen('rho_kgrid_std','w');

f_A = 0.15; f_B = 1 - f_A;
a = 2; c = a*1.633; sigma_smear = 0.53033; %smear R = 0.129099

N_particle = 2; 

%Rsph = (3*f_C*a*a*c*0.43301*4/(4*pi*N_particle))^(1/3);
%Rsph = (3*f_C*a*a*c*0.43301/(4*pi*N_particle))^(1/3);
Rsph = 0.20;

%const = f_C/N_particle/2;
const = f_A/N_particle;
%const = (f_C/N_particle)*2;
%const = 1/(a*a*c*0.866);

form_factor =@(x) 3*(sin(x) - x*cos(x))/x^3;  % x = qR

grid = [30  30  48]; dim = 3;
nodes = (grid(1)/2+1)*grid(2)*grid(3);

rho_A=zeros(nodes,1);
rho_B=zeros(nodes,1);
rho_C=zeros(nodes,1);
wave =zeros(nodes,3);

%file_id = fopen(output_filename,'w');
fprintf(file_id,'  format  1  0\n'); 			% Summary of the system in
fprintf(file_id,'dim\n'); 						% output (rho_kgrid) file
fprintf(file_id,' \t       %d\n',dim);
fprintf(file_id,'crystal_system\n');
fprintf(file_id,' \t       ''%s'' \n','hexagonal');
fprintf(file_id,'N_cell_param\n');
fprintf(file_id,' \t       %i \n',2);
fprintf(file_id,'cell_param\n');
fprintf(file_id,' \t       %6.4f  %6.4f\n',a,c);
fprintf(file_id,'group_name\n');
fprintf(file_id,' \t       ''%s'' \n','p 63/m m c');
fprintf(file_id,'N_monomer\n');
fprintf(file_id,' \t       %d \n',2);
fprintf(file_id,'ngrid\n');
fprintf(file_id,' \t       %i \t\t %i \t\t %i \n',grid(1),grid(2),grid(3));


t = 0;
for i=0:grid(1)/2
    for j=0:grid(2)-1  
        for k=0:grid(3)-1 
           
             G = [i j k];
             G = G_to_bz(G,grid,dim);           % Tranforming waves to first
                                                % brizilion zone (Aliasing)
        
            % Reflection condition for the BCC phase: i+j+k should be even.
            if(G(1)==0 && G(2)==0 && G(3)==0)
               t = t+1;
               wave(t,:)  = G;
               rho_A(t,1) = f_A; 
               rho_B(t,1) = f_B;
              
            else
            % Waves not satisfying the reflection conditions have rho = 0.
               t = t+1;
               wave(t,:) = G;
               
               qR = (((((pi*(2/a))^2)*(G(1))^2)+((((pi*(2/a))^2))*(4/3)*((0.5*G(1))+G(2))^2)+(((pi*(2/c))^2)*G(3)*G(3)))^0.5) * Rsph;
               
               [R,I] = sigma_ff(G(1),G(2),G(3));
               
               rho_A(t,1) = const*form_factor(qR)*exp(-(sigma_smear)^2*qR^2/2)*R;
               rho_B(t,1) = -rho_A(t,1);

            end
            
                           
             fprintf(file_id,'(%6.4E,%6.4E)  (%6.4E,%6.4E)\n', ...
                               rho_A(t,1),0.0,rho_B(t,1),0.0);

       end
    end
end

end

function [R,I] = sigma_ff(h,k,l)

 pos = zeros(1,3);
 
 pos(1,:) = [2/3 1/3 0.75]; 
 pos(2,:) = [1/3 2/3 0.25];


R = 0; I = 0;
for i =1:length(pos(:,1))
    q_dot_r = 2*pi*(h*pos(i,1) + k*pos(i,2) + l*pos(i,3)); 
    R = R + cos(q_dot_r);
    I = I + sin(q_dot_r);
end

end


function Gout = G_to_bz(G,grid,dim)

Gout = G;
if(dim==3)
    if(G(2) > grid(2)/2)
        Gout(2) = G(2) - grid(2);
    end
    
    if(G(3) > grid(3)/2)
        Gout(3) = G(3) - grid(3);
    end
end

if(dim==2)
    if(G(2) > grid(2)/2)
        Gout(2) = G(2) - grid(2);
    end
end

end

