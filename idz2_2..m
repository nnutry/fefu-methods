format short

A = [16  3  1 -2  4;
      2 15 -1  3  2;
      1 -2 17  1  3;
      4  1  2 16 -3;
      3  2  1  4 18];

b = [7; 24; 19; 19; 23];

n = rows(A);

% ---------- Проверка угловых миноров ----------
printf("\n--- Угловые миноры ---\n");
minors_good = true;
for k = 1:n
  d = det(A(1:k, 1:k));
  printf("Delta_%d = %s\n", k, rats(d));
  if d < 1e-12
    minors_good = false;
  endif
endfor

if !minors_good
  error("\nОШИБКА: есть нулевой угловой минор.\n");
endif
printf("Все миноры ненулевые, LU-разложение существует.\n");

% ------- Разложение A = L*U ----------
L = eye(n);
U = zeros(n);

for i = 1:n
  % i-я строка U
  U(i, i:n) = A(i, i:n) - L(i, 1:i-1) * U(1:i-1, i:n);
  % i-й столбец L (ниже диагонали)
  L(i+1:n, i) = (A(i+1:n, i) - L(i+1:n, 1:i-1) * U(1:i-1, i)) / U(i, i);
endfor

printf("\nМатрица L:\n");
disp(L);
printf("\nМатрица U:\n");
disp(U);

printf("\nПроверка: A - L*U должна быть нулевой матрицей:\n");
disp(A - L*U);

% ---------- Ly = b (прямая подстановка) ----------
y = zeros(n, 1);
for i = 1:n
  y(i) = b(i) - L(i, 1:i-1) * y(1:i-1);
endfor
printf("\nВектор y:\n");
disp(y);

% ---------- Ux = y (обратная подстановка) ----------
x = zeros(n, 1);
for i = n:-1:1
  x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
endfor
printf("\nРешение x:\n");
disp(x);

printf("\nНевязка ||Ax - b|| = %e\n", norm(A*x - b));
