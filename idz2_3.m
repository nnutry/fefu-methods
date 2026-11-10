A = [2, 3, 5, -2, 1;
     4, 1, -1, 3, 6;
     1, 5, 2, 4, -2;
     3, -2, 4, 6, 1;
     5, 1, 3, 2, 7];


n = rows(A);
M = [A, eye(n)];     % расширенная матрица

for col = 1:n
  [maxim, idx] = max(abs(M(col:n, col)));
  p = idx + col - 1;

  if maxim < 1e-12
    error("Матрица вырождена, обратной не существует");
  endif

  % перестановка строк
  if p != col
    M([col, p], :) = M([p, col], :);
  endif

  % ведущий элемент делаем единицей
  M(col, :) = M(col, :) / M(col, col);

  % обнуляем столбец во всех остальных строках
  for r = [1:col-1, col+1:n]
    M(r, :) -= M(r, col) * M(col, :);
  endfor
endfor

A_inv = M(:, n+1:2*n);    % правая половина

disp("Обратная матрица:");
disp(A_inv);

disp("A*A_inv:");
disp(A * A_inv);

disp("Погрешность относительно E:");
err = norm(A*A_inv - eye(n), inf)

disp("Сравнение со встроенной inv:");
diff_inv = norm(A_inv - inv(A), inf)
