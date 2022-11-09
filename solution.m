function disp = solution(nDof,fixedDof,K,force)
%fun��o que resolve K*U=F
%fixedDof = Dofs que s�o nulos (n�o precisam ser resolvidos)
%activeDof: Dofs que efetivamente tem valor
%lembrando que disp = displacements (n�o confundir com a fun��o do matlab)
activeDof = setdiff((1:nDof)',fixedDof);
U = K(activeDof,activeDof)\force(activeDof);
disp = zeros(nDof,1);
disp(activeDof) = U
end

