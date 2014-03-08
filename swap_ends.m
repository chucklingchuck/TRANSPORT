function B = swap_ends(A)
  B = A;
  B(:, 1) = A(:, size(A, 2));
  B(:, size(A,2)) = A(:, 1);
end