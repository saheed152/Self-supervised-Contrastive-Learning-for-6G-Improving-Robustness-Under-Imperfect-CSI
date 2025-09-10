clc;
H=H_new(:,:,:,1);
User=5;
desired_user=1;

for i=1:User
    [U,S,V]=svd(H(:,:,i));
    Fopt_i=V(:,[1:2]);
    Wopt_i=U(:,[1:2]);
    Fopt(:,((i-1)*2+1):i*2)=Fopt_i;
    Wopt(:,((i-1)*2+1):i*2)=Wopt_i;
end
[Uk,Sk,Vk]=svd(H(:,:,1));
I=0;
for j=1:User
    if j ~= desired_user
        I=I+Vk'*Fopt(:,[((j-1)*2+1):j*2]);
    end
end

Interference_k_1=Sk([1,2],:)*I;

Interference_k_2=0;
for i=1:User
    if i ~= desired_user
        Interference_k_2= Interference_k_2 + Wopt(:,[1:2])'*H(:,:,1)*Fopt(:,[((i-1)*2+1):i*2]);
    end
end

%% simplified interference 2
Interference_simple=0;
sigma=[Sk(1,1);Sk(2,2)]
for j=1:User
    if j ~= desired_user
        Interference_simple=Interference_simple+sigma.*Vk(:,[1,2])'*Fopt(:,[((j-1)*2+1):j*2]);
    end
end


%% comparison between previous intefernce calculation and the actual one
interference_before= sigma.*Vk(:,[1,2])'*Fopt(:,[3:10]);

%% More simplification of (2): 4
for i=1:User
    [U,S,V]=svd(H(:,:,i));
    Fopt_i=V(:,[1:2]);
    Fopt2(((i-1)*144+1):i*144,:)=Fopt_i;
end
for i=1:User
    [U,S,V]=svd(H(:,:,1));
    if i ~= desired_user
        Fopt_i=V(:,[1:2])';
        Vp(:,((i-1)*144+1):i*144)=Fopt_i;

    end
end

interference_after = sigma.*Vp*Fopt2;