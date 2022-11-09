function disp = solution(nDof,fixedDof,K,force)
%função que resolve K*U=F
%fixedDof = Dofs que são nulos (não precisam ser resolvidos)
%activeDof: Dofs que efetivamente tem valor
%lembrando que disp = displacements (não confundir com a função do matlab)
activeDof = setdiff((1:nDof)',fixedDof);
U = K(activeDof,activeDof)\force(activeDof);
disp = zeros(nDof,1);
disp(activeDof) = U
end

