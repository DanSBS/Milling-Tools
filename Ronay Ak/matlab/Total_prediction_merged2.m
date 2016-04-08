clear all

tic
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,fliplr(lower)],color),'EdgeColor',color);
process_ID=[1 2 3 4 5 6]; %[cut_with_feed, plung, air-cut, rapid-motion, dwell]

%base covariance function scale parameters
Hyp = [log(1000),log(6000),log(6),log(3),log(3),log(3),log(3),log(3),log(3),log(3),log(1)];

     
Feature_set{1} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{2} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{3} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{4} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{4} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{5} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{6} = [1 2 3   4 5 6 7   8 9 10 ];




load('Training1.mat')
load('Training2.mat')
load('Training3.mat')
load('Training4.mat')
load('Training5.mat')
load('Training6.mat')
load('Training7.mat')
load('Training8.mat')
load('Training9.mat')
load('Training10.mat')
load('Training11.mat')
load('Training12.mat')
load('Training13.mat')
load('Training14.mat')
load('Training15.mat')
load('Training16.mat')
load('Training17.mat')
load('Training18.mat')

D = [Training1;Training2;Training3;Training4;Training5;Training6;Training7;Training8;Training9;Training10;Training11;Training12;Training13;Training14;Training15;Training16;Training17;Training18];
%D = [Training1;Training2;Training3;Training4;Training5;Training6;Training7;Training8;Training9];
%D = [Training10;Training11;Training12;Training13;Training14;Training15;Training16;Training17;Training18];

%% extract the field
energy = cell2mat(D(:,9)); %energy consumption
duration = cell2mat(D(:,10)); %duration of operation
feed = cell2mat(D(:,11)); %duration of operation
spindle = cell2mat(D(:,12)); %spindle speed
length_cut_X = abs(cell2mat(D(:,19))); %code dx
length_cut_Y = abs(cell2mat(D(:,20))); %code dy
length_cut_Z = abs(cell2mat(D(:,23))); %code dy
length_cut_XY = cell2mat(D(:,21)); %code length_cut
length_cut_XYZ = cell2mat(D(:,31)); %code length_cut
actual_dx = cell2mat(D(:,25)); %actual dx 
actual_dy = cell2mat(D(:,26)); %actual dy
actual_length_cut = sqrt(actual_dx.^2+actual_dy.^2);
depth_cut = cell2mat(D(:,27)); %Depth of cut
area_cut = cell2mat(D(:,28)); %Depth of cut
volume_cut = cell2mat(D(:,29)); %Depth of cut

%ratio cut
ratio_cut = actual_length_cut./length_cut_XY;
for i=1:length(ratio_cut)
    if ratio_cut(i)>0;
        ratio_cut(i)=ratio_cut(i);
    else
        ratio_cut(i)=0;
    end
end

%cut_method=zeros(length(y),1);
for i=1:length(energy)

    if strcmp(D{i,24},'Conventional')
        cut_method(i,:) = [1,0,0];   
    elseif strcmp(D{i,24},'Climb')
        cut_method(i,:) = [0,1,0];
    else
        cut_method(i,:)=[0,0,1];
    end 
end

%cut_direction=zeros(length(y),1);
for i=1:length(energy)
    if abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))==0
        cut_direction(i,:) = [1,0,0,0];   
    elseif abs(cell2mat(D(i,19)))==0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,:) = [0,1,0,0]; 
    elseif abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,:)=[0,0,1,0]; 
    else
        cut_direction(i,:)=[0,0,0,1]; 
    end
        
end





%label
for i=1:length(energy)
    if strcmp(D{i,30},'Cut with Feed')
        label(i,1) = 1;
    elseif strcmp(D{i,30},'Plunge with feed')
        label(i,1) = 2;
        
        
        
    elseif strcmp(D{i,30},'Air-Cut')
        label(i,1) = 3;
    elseif strcmp(D{i,30},'Air-Cut in Z while plunging')
        label(i,1) = 4;
    elseif strcmp(D{i,30},'Air-cut in Z while retracting')
        label(i,1) = 5;
        
        
        
    elseif strcmp(D{i,30},'No Cut - Rapid motion')
        label(i,1) = 6;
           
    elseif strcmp(D{i,30},'Dwell')
        label(i,1) = 7;  
        
    else
        label(i,1) = 0;
    end   
end




    for i=1:length(energy)
        % cut with feed
        if ( label(i)==1 ) %cut with feed (x, y z included)
            ID(i,1) = 1;
        elseif ( label(i)==2)% plunge with feed
            ID(i,1) = 2;
        elseif ( label(i)==3 || label(i)==4 || label(i)==5)% air cut in x-y
            ID(i,1) = 3;
        elseif (label(i)==6) %Rapid motion
            ID(i,1)= 4; 
        elseif (label(i)==7) %dwell
            ID(i,1)= 5; 
        else
            ID(i,1) = 6;
        end
    end
    
    
input= [feed,...,        %1
     spindle,...,        %2
     depth_cut,...,      %3
     cut_direction,...,  %4 5 6 7
     cut_method,...,     %8 9 10
     length_cut_XYZ,..., %11
     ID,...,             %12
     duration];          %13
output =energy;
%density=energy./length_cut_XYZ;

%total data set
X = input;
E = energy;
for i=1:length(E)
    if ID(i) == 5 %dwell
        L(i,1)=1;
    else
        L(i,1)=X(i,11);  
    end
end
Y = E./L;

index_zero = find(Y==0);
index_zero=union([index_zero], [index_zero+1])
index_zero=union([index_zero], [index_zero-1])

%clean up data
clean_up_index = setdiff([1:length(energy)],[index_zero]);
X = X(clean_up_index,:);
ID = ID(clean_up_index);
L = L(clean_up_index);
E = E(clean_up_index);
Y = Y(clean_up_index);



%clean up data
clean_up_index = find(Y > 0 & Y < 3000);
X = X(clean_up_index,:);
ID = ID(clean_up_index);
L = L(clean_up_index);
E = E(clean_up_index);
Y = Y(clean_up_index);




%training individul prediction function
for I=1:length(process_ID)
    
    feature_index = Feature_set{I}
    process_id = process_ID(I);
    
    %Select data depending on the process ID
    if process_id == 4 %rapid motion
        index_job{I} = find(X(:,12)==process_id );
    else                                             %time           %non-zero cut
        index_job{I} = find(X(:,12)==process_id & X(:,13)>2 & X(:,11)>0); %cutting related (1 is the best)
    end

    %selection data
    X_training=X(index_job{I},feature_index);
    E_training=E(index_job{I});
    L_training=L(index_job{I});
    T_training=X(index_job{I},13);
    Y_training=E_training./L_training;

    %Learning GP
    hyp.cov = [Hyp(feature_index),0];
    meanfunc = {@meanSum, {@meanLinear, @meanConst}}; %hyp.mean = zeros(size(X,2)+1,1);
    likfunc = @likGauss; sn = 500; hyp.lik = log(sn);
    covfunc = @covSEard; 
    hyp2 = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, X_training, Y_training)
    f{I}=hyp2
    
end


toc










inputFiles_blind{1}='blind_test_accurate.mat'
inputFiles_blind{2}='blind_test_intermediate.mat'
inputFiles_blind{3}='blind_test_bad.mat'


for I=1:length(inputFiles_blind)

    %prediction part
    clear cut_method cut_direction type_operation label ID input E_predict S_predict Y_predict L_blind Y_blind X_blind
    load(inputFiles_blind{I})
    
    D = [Data];

    
    
    %extract fetures
    energy = cell2mat(D(:,9)); %answer
    duration = cell2mat(D(:,10)); %duration
    feed = cell2mat(D(:,11)); %duration of operation
    spindle = cell2mat(D(:,12)); %spindle speed
    length_cut_X = abs(cell2mat(D(:,19))); %code dx
    length_cut_Y = abs(cell2mat(D(:,20))); %code dy
    length_cut_Z = abs(cell2mat(D(:,23))); %code dy
    length_cut_XY = cell2mat(D(:,21)); %code length_cut
    length_cut_XYZ = cell2mat(D(:,31)); %code length_cut
    actual_dx = cell2mat(D(:,25)); %actual dx 
    actual_dy = cell2mat(D(:,26)); %actual dy
    actual_length_cut = sqrt(actual_dx.^2+actual_dy.^2);
    depth_cut = cell2mat(D(:,27)); %Depth of cut
    volume_cut = cell2mat(D(:,29)); %Depth of cut

    ratio_cut = actual_length_cut./length_cut_XY;
    for i=1:length(ratio_cut)
        if ratio_cut(i)>0;
            ratio_cut(i)=ratio_cut(i);
        else
            ratio_cut(i)=0;
        end
    end

    %cut_method=zeros(length(y),1);
for i=1:length(energy)

    if strcmp(D{i,24},'Conventional')
        cut_method(i,:) = [1,0,0];   
    elseif strcmp(D{i,24},'Climb')
        cut_method(i,:) = [0,1,0];
    else
        cut_method(i,:)=[0,0,1];
    end 
end

%cut_direction=zeros(length(y),1);
for i=1:length(energy)
    if abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))==0
        cut_direction(i,:) = [1,0,0,0];   
    elseif abs(cell2mat(D(i,19)))==0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,:) = [0,1,0,0]; 
    elseif abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,:)=[0,0,1,0]; 
    else
        cut_direction(i,:)=[0,0,0,1]; 
    end
        
end


    %label
    for i=1:length(feed)

        if strcmp(D{i,30},'Cut with Feed')
            label(i,1) = 1;
        elseif strcmp(D{i,30},'Plunge with feed')
            label(i,1) = 2;



        elseif strcmp(D{i,30},'Air-Cut')
            label(i,1) = 3;
        elseif strcmp(D{i,30},'Air-Cut in Z while plunging')
            label(i,1) = 4;
        elseif strcmp(D{i,30},'Air-cut in Z while retracting')
            label(i,1) = 5;


        elseif strcmp(D{i,30},'No Cut - Rapid motion')
            label(i,1) = 6;
        elseif strcmp(D{i,30},'Rapid retract')
            label(i,1) = 6;         

        elseif strcmp(D{i,30},'Dwell')
            label(i,1) = 7;  
        else
            label(i,1) = 0;
        end   
    end


    for i=1:length(energy)
        % cut with feed
        if ( label(i)==1 ) %cut with feed (x, y z included)
            ID(i,1) = 1;
        elseif ( label(i)==2)% plunge with feed
            ID(i,1) = 2;
        elseif ( label(i)==3 || label(i)==4 || label(i)==5)% air cut in x-y
            ID(i,1) = 3;
        elseif (label(i)==6) %Rapid motion
            ID(i,1)= 4; 
        elseif (label(i)==7) %dwell
            ID(i,1)= 5; 
        else
            ID(i,1) = 6;
        end
    end



input= [feed,...,        %1
     spindle,...,        %2
     depth_cut,...,      %3
     cut_direction,...,  %4 5 6 7
     cut_method,...,     %8 9 10
     length_cut_XYZ,..., %11
     ID,...,             %12
     duration];          %13

    output =energy;
    density=energy./length_cut_XYZ;


    %total data set
    X_blind = input;
    ID_blind = input(:,12);
    E_blind = energy;
    for i=1:length(E_blind)
        if ID(i) == 5 %dwell
            L_blind(i,1)=1;
        else
            L_blind(i,1)=X_blind(i,11);  
        end
    end
    Y_blind = E_blind./L_blind;
    T_blind = input(:,13);

    
    
    index_zero = find(Y_blind==0);
    index_zero=union([index_zero], [index_zero+1])
    index_zero=union([index_zero], [index_zero-1])

    %clean up data
    clean_up_index = setdiff([1:length(energy)],[index_zero]);
    X_blind = X_blind(clean_up_index,:);
    ID_blind = ID_blind(clean_up_index);
    L_blind = L_blind(clean_up_index);
    E_blind = E_blind(clean_up_index);
    Y_blind = Y_blind(clean_up_index);
    T_blind = T_blind(clean_up_index);  
    

    % %clean up data
    clean_up_index = find(Y_blind > 0 & Y_blind < inf);
    X_blind = X_blind(clean_up_index,:);
    ID_blind = ID_blind(clean_up_index);
    L_blind = L_blind(clean_up_index);
    E_blind = E_blind(clean_up_index);
    Y_blind = Y_blind(clean_up_index);
    T_blind = T_blind(clean_up_index);





    for i=1:length(E_blind)


            %select the right prediction function

            process_id = ID_blind(i); 

            if (process_id>0 && process_id <6)
                hyp = f{process_id};
                X_training = X(index_job{process_id}, Feature_set{process_id});
                Y_training = Y(index_job{process_id});

                [Y_predict(i,1), S_predict(i,1) ] = gp(hyp, @infExact, [], covfunc, likfunc, X_training, Y_training, X_blind(i,Feature_set{process_id}));
                %S_predict(i,1)=S_predict(i,1)+sn;
                E_predict(i,1) = Y_predict(i,1).*L_blind(i);
            else
                Y_predict(i,1)=0;
                S_predict(i,1)=0;
                E_predict(i,1)=0;     
%                 Y_blind(i,1)=0;
%                 S_blind(i,1)=0;
%                 E_blind(i,1)=0;               

            end
    end

    Number_b{I} = length(T_blind);
    Duration_b{I} = mean(T_blind);
    RAE_density_b{I} = mean(abs(Y_predict.*L_blind-Y_blind.*L_blind))/mean(Y_blind.*L_blind);
    predicted_energy_b{I} = sum(Y_predict.*L_blind);
    measured_energy_b{I} = sum(E_blind);
    standard_b{I} = sqrt(sum(S_predict.*L_blind.^2));
    RTE_density_b{I} = (sum(Y_predict.*L_blind)-sum(E_blind))/sum(E_blind);

    %display([number,duration,RAE_density*100,predicted_energy/1000,measured_energy/1000,standard,RTE_density*100])

    YY_measured_b{I} = Y_blind;
    EE_measured_b{I} = E_blind;
    LL_measured_b{I} = L_blind;

    EE_predict_b{I} = E_predict;
    SS_predict_b{I} = S_predict;
    YY_predict_b{I} = Y_predict;





end


display([Number_b{1},Duration_b{1},RAE_density_b{1}*100,measured_energy_b{1}/1000, predicted_energy_b{1}/1000,standard_b{1}/1000,RTE_density_b{1}*100;
        Number_b{2},Duration_b{2},RAE_density_b{2}*100,measured_energy_b{2}/1000,predicted_energy_b{2}/1000,standard_b{2}/1000,RTE_density_b{2}*100;
        Number_b{3},Duration_b{3},RAE_density_b{3}*100,measured_energy_b{3}/1000,predicted_energy_b{3}/1000,standard_b{3}/1000,RTE_density_b{3}*100])


figure(1)
subplot(3,1,1)
hold
stairs(EE_measured_b{1},'r')
stairs(EE_predict_b{1},'b')
ylabel('Energy (J)','Interpreter','Latex')
subplot(3,1,2)
hold
stairs(EE_measured_b{2},'r')
stairs(EE_predict_b{2},'b')
ylabel('Energy (J)','Interpreter','Latex')
subplot(3,1,3)
hold
stairs(EE_measured_b{3},'r')
stairs(EE_predict_b{3},'b')
ylabel('Energy (J)','Interpreter','Latex')





























inputFiles{1}='toolPath1.mat'
inputFiles{2}='toolPath2.mat'
inputFiles{3}='toolPath3.mat'
inputFiles{4}='toolPath4.mat'

for I=1:length(inputFiles)

    %prediction part
    clear cut_method cut_direction type_operation label ID input E_predict S_predict Y_predict L_blind Y_blind X_blind
    load(inputFiles{I})
    if I==1
        D = [toolPath1];
    elseif I==2
        D = [toolPath2];
    elseif I==3
        D = [toolPath3];
    else
        D = [toolPath4];
    end

    
    
    %extract fetures
    energy = cell2mat(D(:,9)); %answer
    duration = cell2mat(D(:,10)); %duration
    feed = cell2mat(D(:,11)); %duration of operation
    spindle = cell2mat(D(:,12)); %spindle speed
    length_cut_X = abs(cell2mat(D(:,19))); %code dx
    length_cut_Y = abs(cell2mat(D(:,20))); %code dy
    length_cut_Z = abs(cell2mat(D(:,23))); %code dy
    length_cut_XY = cell2mat(D(:,21)); %code length_cut
    length_cut_XYZ = cell2mat(D(:,31)); %code length_cut
    actual_dx = cell2mat(D(:,25)); %actual dx 
    actual_dy = cell2mat(D(:,26)); %actual dy
    actual_length_cut = sqrt(actual_dx.^2+actual_dy.^2);
    depth_cut = cell2mat(D(:,27)); %Depth of cut
    volume_cut = cell2mat(D(:,29)); %Depth of cut

    ratio_cut = actual_length_cut./length_cut_XY;
    for i=1:length(ratio_cut)
        if ratio_cut(i)>0;
            ratio_cut(i)=ratio_cut(i);
        else
            ratio_cut(i)=0;
        end
    end

    %cut_method=zeros(length(y),1);
for i=1:length(energy)

    if strcmp(D{i,24},'Conventional')
        cut_method(i,:) = [1,0,0];   
    elseif strcmp(D{i,24},'Climb')
        cut_method(i,:) = [0,1,0];
    else
        cut_method(i,:)=[0,0,1];
    end 
end

%cut_direction=zeros(length(y),1);
for i=1:length(energy)
    if abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))==0
        cut_direction(i,:) = [1,0,0,0];   
    elseif abs(cell2mat(D(i,19)))==0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,:) = [0,1,0,0]; 
    elseif abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,:)=[0,0,1,0]; 
    else
        cut_direction(i,:)=[0,0,0,1]; 
    end
        
end


    %opeartion 
    for i=1:length(feed)

        if strcmp(D{i,33},'Face Milling')
            type_operation(i,1) = 1;      
        elseif strcmp(D{i,33},'Contouring')
            type_operation(i,1) = 2;       
        elseif strcmp(D{i,33},'Slotting')
            type_operation(i,1) = 3;
        elseif strcmp(D{i,33},'Pocketing')
            type_operation(i,1) = 4;
        elseif strcmp(D{i,33},'Spiraling')
            type_operation(i,1) = 5;
        elseif strcmp(D{i,33},'Drilling')
            type_operation(i,1) = 6;             
        else
            type_operation(i,1) = 0;
        end    
    end



    %label
    for i=1:length(feed)

        if strcmp(D{i,30},'Cut with Feed')
            label(i,1) = 1;
        elseif strcmp(D{i,30},'Plunge with feed')
            label(i,1) = 2;



        elseif strcmp(D{i,30},'Air-Cut')
            label(i,1) = 3;
        elseif strcmp(D{i,30},'Air-Cut in Z while plunging')
            label(i,1) = 4;
        elseif strcmp(D{i,30},'Air-cut in Z while retracting')
            label(i,1) = 5;


        elseif strcmp(D{i,30},'No Cut - Rapid motion')
            label(i,1) = 6;
        elseif strcmp(D{i,30},'Rapid retract')
            label(i,1) = 6;         

        elseif strcmp(D{i,30},'Dwell')
            label(i,1) = 7;  
        else
            label(i,1) = 0;
        end   
    end

    for i=1:length(energy)
        % cut with feed
        if ( label(i)==1 ) %cut with feed (x, y z included)
            ID(i,1) = 1;
        elseif ( label(i)==2)% plunge with feed
            ID(i,1) = 2;
        elseif ( label(i)==3 || label(i)==4 || label(i)==5)% air cut in x-y
            ID(i,1) = 3;
%         elseif (label(i)==5 ) %air cut in z
%             ID(i,1)= 2; 
        elseif (label(i)==6) %Rapid motion
            ID(i,1)= 4; 
        elseif (label(i)==7) %dwell
            ID(i,1)= 5; 
        else
            ID(i,1) = 6;
        end
    end


input= [feed,...,        %1
     spindle,...,        %2
     depth_cut,...,      %3
     cut_direction,...,  %4 5 6 7
     cut_method,...,     %8 9 10
     length_cut_XYZ,..., %11
     ID,...,             %12
     duration];  
 
    output =energy;
    density=energy./length_cut_XYZ;


    %total data set
    X_blind = input;
    ID_blind = input(:,12);
    E_blind = energy;
    for i=1:length(E_blind)
        if ID(i) == 5 %dwell
            L_blind(i,1)=1;
        else
            L_blind(i,1)=X_blind(i,11);  
        end
    end
    Y_blind = E_blind./L_blind;
    T_blind = input(:,13);


    index_zero = find(Y_blind==0);
    index_zero=union([index_zero], [index_zero+1])
    index_zero=union([index_zero], [index_zero-1])

    %clean up data
    clean_up_index = setdiff([1:length(energy)],[index_zero]);
    X_blind = X_blind(clean_up_index,:);
    ID_blind = ID_blind(clean_up_index);
    L_blind = L_blind(clean_up_index);
    E_blind = E_blind(clean_up_index);
    Y_blind = Y_blind(clean_up_index);
    T_blind = T_blind(clean_up_index);      
    
    
    % %clean up data
    clean_up_index = find(Y_blind > 0 & Y_blind < inf);
    X_blind = X_blind(clean_up_index,:);
    ID_blind = ID_blind(clean_up_index);
    L_blind = L_blind(clean_up_index);
    E_blind = E_blind(clean_up_index);
    Y_blind = Y_blind(clean_up_index);
    T_blind = T_blind(clean_up_index);





    for i=1:length(E_blind)


            %select the right prediction function

            process_id = ID_blind(i); 

            if (process_id>0 && process_id <6)
                hyp = f{process_id};
                X_training = X(index_job{process_id}, Feature_set{process_id});
                Y_training = Y(index_job{process_id});

                [ Y_predict(i,1), S_predict(i,1) ] = gp(hyp, @infExact, [], covfunc, likfunc, X_training, Y_training, X_blind(i,Feature_set{process_id}));
                %S_predict(i,1)=S_predict(i,1)+sn;
                E_predict(i,1) = Y_predict(i,1).*L_blind(i);
            else
                %Y_predict(i,1)=0;
                %S_predict(i,1)=0;
                %E_predict(i,1)=0;     
%                 Y_blind(i,1)=0;
%                 S_blind(i,1)=0;
%                 E_blind(i,1)=0;               

            end
    end

    Number_t{I} = length(T_blind);
    Duration_t{I} = mean(T_blind);
    RAE_density_t{I} = mean(abs(Y_predict.*L_blind-Y_blind.*L_blind))/mean(Y_blind.*L_blind);
    predicted_energy_t{I} = sum(Y_predict.*L_blind);
    measured_energy_t{I} = sum(E_blind);
    standard_t{I} = sqrt(sum(S_predict.*L_blind.^2));
    RTE_density_t{I} = (sum(Y_predict.*L_blind)-sum(E_blind))/sum(E_blind);

    %display([number,duration,RAE_density*100,predicted_energy/1000,measured_energy/1000,standard,RTE_density*100])
    XX_measured_t{I} = X_blind;
    YY_measured_t{I} = Y_blind;
    EE_measured_t{I} = E_blind;
    LL_measured_t{I} = L_blind;

    EE_predict_t{I} = E_predict;
    SS_predict_t{I} = S_predict;
    YY_predict_t{I} = Y_predict;





end



figure(2)
subplot(4,1,1)
hold
stairs(EE_measured_t{1},'r')
stairs(EE_predict_t{1},'b')
ylabel('Energy (J)','Interpreter','Latex')
subplot(4,1,2)
hold
stairs(EE_measured_t{2},'r')
stairs(EE_predict_t{2},'b')
ylabel('Energy (J)','Interpreter','Latex')
subplot(4,1,3)
hold
stairs(EE_measured_t{3},'r')
stairs(EE_predict_t{3},'b')
ylabel('Energy (J)','Interpreter','Latex')
subplot(4,1,4)
hold
stairs(EE_measured_t{4},'r')
stairs(EE_predict_t{4},'b')
ylabel('Energy (J)','Interpreter','Latex')



display([Number_t{1},Duration_t{1},RAE_density_t{1}*100,measured_energy_t{1}/1000,predicted_energy_t{1}/1000,standard_t{1}/1000,RTE_density_t{1}*100;
        Number_t{2},Duration_t{2},RAE_density_t{2}*100,measured_energy_t{2}/1000,predicted_energy_t{2}/1000,standard_t{2}/1000,RTE_density_t{2}*100;
        Number_t{3},Duration_t{3},RAE_density_t{3}*100,measured_energy_t{3}/1000,predicted_energy_t{3}/1000,standard_t{3}/1000,RTE_density_t{3}*100;
        Number_t{4},Duration_t{4},RAE_density_t{4}*100,measured_energy_t{4}/1000,predicted_energy_t{4}/1000,standard_t{4}/1000,RTE_density_t{4}*100])




figure(3)
clear label
y = [sum(EE_measured_t{1}) sum(EE_predict_t{1});
     sum(EE_measured_t{2}) sum(EE_predict_t{2});
     sum(EE_measured_t{3}) sum(EE_predict_t{3});
     sum(EE_measured_t{4}) sum(EE_predict_t{4})]/1000;

bound = [standard_t{1};standard_t{2};standard_t{3};standard_t{4}]/1000;

errY=zeros(4,2,2);
errY(:,2,1)=bound*1.96;
errY(:,2,2)=bound*1.96;
barwitherr(errY, y) 

%bar(y)
set(gca,'XTickLabel',{'Path 1','Path 2','Path 3', 'Path 4'})
ylabel('Energy (kJ)','Interpreter','Latex')
legend('Measured','Predicted')
ylim([0,420])




load('toolPaths.mat')
figure(4)
subplot(1,4,1)
hold
plot(P1(4,1),P1(4,2),'ro','Markerfacecolor','r','Markersize',8)
legend('Start')
plot(P1(4:end,1),P1(4:end,2))
xlabel('x-axis')
ylabel('y-axis')

subplot(1,4,2)
hold
plot(P2(4,1),P2(4,2),'ro','Markerfacecolor','r','Markersize',8)
legend('Start')
plot(P2(4:end,1),P2(4:end,2))
xlabel('x-axis')
ylabel('y-axis')

subplot(1,4,3)
hold
plot(P3(4,1),P3(4,2),'ro','Markerfacecolor','r','Markersize',8)
legend('Start')
plot(P3(4:end,1),P3(4:end,2))
xlabel('x-axis')
ylabel('y-axis')

subplot(1,4,4)
hold
plot(P4(4,1),P4(4,2),'ro','Markerfacecolor','r','Markersize',8)
legend('Start')
plot(P4(4:end,1),P4(4:end,2))
xlabel('x-axis')
ylabel('y-axis')


