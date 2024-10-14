%Compatibility Calculator - Calculates the compatability of two individuals
%based on their fields and levels of interest.
clear all;
t = 0:(pi/24):12;

% Gather some information
%Person 1
fprintf('For the following categories, please select your prefered taste, \nas well as a level of importance\n');
p1 = input('Movies:\n1. Horror\n2. Comedy\n3. Action\n4. Drama\n5. Documentary\n');
a1 = input('\nEnter a value from 1 - 10 for your interest level:\n');
p2 = input('\nMusic:\n1. Classical\n2. Blues\n3. Rock\n4. Pop\n5. Rap\n');
a2 = input('\nEnter a value from 1 - 10 for your interest level:\n');
p3 = input('\nSports: \n1. Soccer\n2. Football\n3. Baseball\n4. Cheerleading\n5. Basketball\n');
a3 = input('\nEnter a value from 1 - 10 for your interest level:\n');
p4 = input('\nReligious Beliefs: \n1. Christianity \n2. Judaism \n3. Islam \n4. Spiritual \n5. Not Religiouslly Affiliated\n');
a4 = input('\nEnter a value from 1 - 10 for your involvment level:\n');
p5 = input('\nExercise Types: \n1. Hiking \n2. Swimming \n3. Jogging \n4. Yoga \n5. None\n');
a5 = input('\nEnter a value from 0 - 7 for how many days a week you exercise:\n');
if (a5 == 0)
    a5 = 9000;
end
p6 = input('\nNews Interests: \n1. Art \n2. Politics \n3. Sports \n4. Technology \n5. Pop Culture\n');
a6 = input('\nEnter a value from 1 - 10 for your interest level:\n');
a7 = input('\nPolitical Affiliation: \n1. Libertarian \n2. Republican \n3. Democrat \n4. Socialist \n5. Other\n');
p7 = input('\nEnter a value from 1 10 for your political involvment: \n');


fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');

%Person 2
fprintf('For the following areas of interest, please select your prefered taste, \nas well as a level of importance\n');
p21 = input('\nMovies:\n1. Horror\n2. Comedy\n3. Action\n4. Drama\n5. Documentary\n');
a21 = input('\nEnter a value from 1 - 10 for your interest level:\n');
p22 = input('\nMusic:\n1. Classical\n2. Blues\n3. Rock\n4. Pop\n5. Rap\n');
a22 = input('\nEnter a value from 1 - 10 for your interest level:\n');
p23 = input('\nSports: \n1. Soccer\n2. Football\n3. Baseball\n4. Cheerleading\n5. Basketball\n');
a23 = input('\nEnter a value from 1 - 10 for your interest level:\n');
p24 = input('\nReligious Beliefs: \n1. Christianity \n2. Judaism \n3. Islam \n4. Spiritual \n5. Not Religiouslly Affiliated\n');
a24 = input('\nEnter a value from 1 - 10 for your involvment level:\n');
p25 = input('\nExercise Types: \n1. Hiking \n2. Swimming \n3. Jogging \n4. Yoga \n5. None\n');
a25 = input('\nEnter a value from 0 - 7 for how many days a week you exercise:\n');
if (a25 == 0)
    a25 = 9000;
end
p26 = input('\nNews Interests: \n1. Art \n2. Politics \n3. Sports \n4. Technology \n5. Pop Culture\n');
a26 = input('\nEnter a value from 1 - 10 for your interest level:\n');
a27 = input('\nPolitical Affiliation: \n1. Libertarian \n2. Republican \n3. Democrat \n4. Socialist \n5. Other\n');
p27 = input('\nEnter a value from 1 10 for your political involvment: \n');

%Wave Calculations
%Person 1
wave1 = sin(t);
wave2 = (1/a1)*sin(2*t*p1);
wave3 = (1/a2)*sin(3*t*p2);
wave4 = (1/a3)*sin(4*t*p3);
wave5 = (1/a4)*sin(5*t*p4);
wave6 = (1/a5)*sin(6*t*p5);
wave7 = (1/a6)*sin(7*t*p6);
wave8 = (1/a7)*sin(8*t*p7);
wave_total1 = wave1 + wave2 + wave3 + wave4 + wave5 + wave6 + wave7 + wave8;

%Person 2
wave21 = sin(t);
wave22 = (1/a21)*sin(2*t*p21);
wave23 = (1/a22)*sin(3*t*p22);
wave24 = (1/a23)*sin(4*t*p23);
wave25 = (1/a24)*sin(5*t*p24);
wave26 = (1/a25)*sin(6*t*p25);
wave27 = (1/a26)*sin(7*t*p26);
wave28 = (1/a27)*sin(8*t*p27);
wave_total2 = wave21 + wave22 + wave23 + wave24 + wave25 + wave26 + wave27 + wave28;

wave_total3 = wave_total2 + wave_total1;

%Graph Outputs
%Person 1
figure(1); clf;
plot(t, wave_total1, 'red');

%Person 2
figure(2); clf;
plot(t, wave_total2, 'blue');

%Total
figure(3); clf;
plot(t, wave_total3, 'green');












