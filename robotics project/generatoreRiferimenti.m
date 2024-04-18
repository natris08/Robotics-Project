L = [1.0; 0.8];
q1 = 0;
q2 = 0;
P1 = [0.5; 0.3];
P2 = [1.5; 0.3];
N = 100;
P1P2 = zeros(2, N);
qq1 = zeros(1, 300);
qq2 = zeros(1, 300);
for i=1:N
    [x, y] = leggeOrariaRetta(0, 10, P1, P2, i);
    P1P2(:,i) = [x, y]';
    Q_attuale = twolink_inverse(L, x, y);
    q1 = Q_attuale(1,1);
    qq1(i) = Q_attuale(1,1);
    qq1(i+200) = Q_attuale(1,1);
    q2 = Q_attuale(2,1);
    qq2(i) = Q_attuale(2,1);
    qq2(i+200) = Q_attuale(2,1);
end

P2P1 = zeros(2, N);
for i=1:N
    Q1_iniziale = twolink_inverse(L, P2(1), P2(2));
    Q2_finale = twolink_inverse(L, P1(1), P1(2));
    [q1, q2] = leggeOrariaSemicirconferenza(0, 10, Q1_iniziale, Q2_finale, i);
    Q = [q1, q2]';
    qq1(i+100) = Q(1);
    qq2(i+100) = Q(2);
    P2P1(:,i) = twolink_direct(L, Q);
end
Q = [qq1;qq2];
tempo = 0.1:0.1:30;
figure(1)
plot(tempo, Q)
legend('q1','q2')
ylabel('rad')
xlabel('Tempo')
QV = gradient(Q);
figure(2)
plot(tempo, QV)
legend('q1','q2')
ylabel('rad')
xlabel('Tempo')
figure(3)
plot(P1P2(1,:), P1P2(2,:),'ro', P2P1(1,:), P2P1(2,:), 'b>', P1P2(1,:), P1P2(2,:),'b>')
legend('P1->P2');

function P=twolink_direct(L, Q)
L1 = L(1);
L2 = L(2);
q1 = Q(1);
q2 = Q(2);
Px = L1*cos(q1)+L2*cos(q1+q2);
Py = L1*sin(q1)+L2*sin(q1+q2);
P = [Px; Py];
end

function Q = twolink_inverse(L, x, y)
L1 = L(1);
L2 = L(2);
Px = x;
Py = y;
P = [Px; Py];
c2 = (Px^2+Py^2-L1^2-L2^2)/(2*L1*L2);
s2 = sqrt(1-c2^2);
z = (1/(L1^2+L2^2+2*L1*L2*c2))*[L1+L2*c2 L2*s2; -L2*s2 L1+L2*c2]*P;
c1 = z(1, 1);
s1 = z(2, 1);
q1 = atan2(s1, c1);
q2 = atan2(s2, c2);
Q = [q1; q2];
end

function [px, py] = leggeOrariaRetta(t1, t2, P1, P2, i)
lambda = calcolaLambda(calcolaSigma(t1, t2, i));
px = P1(1) + lambda*(P2(1) - P1(1));
py = P1(2) + lambda*(P2(2) - P1(2));
end

function [q1, q2] = leggeOrariaSemicirconferenza(t1, t2, Q1, Q2, i)
lambda = calcolaLambda(calcolaSigma(t1, t2, i));
q1 = Q1(1) + lambda * (Q2(1) - Q1(1));
q2 = Q1(2) + lambda * (Q2(2) - Q1(2));
end

function [sigma] = calcolaSigma(t1, t2, i)
t = t1:0.1:t2;
sigma = ((t(i))-t1)/(t2-t1);
end

function [lambda] = calcolaLambda(sigma)
% lambda = -2*sigma^3 + 3*sigma^2;
% lambda = 6*sigma^5 - 15 *sigma^4 + 10*sigma^3;
lambda = -20*sigma^7 + 70*sigma^6 - 84*sigma^5 + 35*sigma^4;
end