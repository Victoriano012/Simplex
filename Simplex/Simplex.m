function [] = simplex()
    fileID = fopen('Dades/PL4_2.txt','r');
    formatSpec = '%d %f';
    n = fscanf(fileID,formatSpec,1).';
    m = fscanf(fileID,formatSpec,1).';
    sizec = [1 n];
    sizeA = [n m];
    sizeb = [1 m];
    c = fscanf(fileID,formatSpec,sizec).';
    A = fscanf(fileID,formatSpec,sizeA).';
    b = fscanf(fileID,formatSpec,sizeb).';
    
    Bland = false;
    
    if Bland
        fprintf('Inici simplex primal amb regla de Bland \n');
    else
        fprintf('Inici simplex primal amb la regla del cost reduït més negatiu \n');
    end
    
    
    % aquí comença l'algorisme del símplex

    % declaració de les variables de la faseI:
    A_I = [A, eye(m)];
    c_I = [zeros(n, 1); ones(m, 1)];
    xb_I = b;
    VB_I = n+1:m+n;
    VNB_I = 1:n;
    invB_I = eye(m);
    
    %faseI del simplex:
    fprintf('  Fase I \n');
    [xb, ~, VB, VNB, invB, ite] = fase(A_I, c_I, xb_I, VB_I, VNB_I, invB_I, 1, Bland);
    fprintf('  q = %2d,   rq = % 9.3f,   B(p) = %2d,   theta* = % 9.3f,   z = % 9.3f \n', 0, 0, 0, 0, c_I(VB)'*xb);
    
    % quan la solució òptima de la faseI sigui diferent de 0 sabem que tenim un problema infactible
    z = c_I(VB)'*xb;
    if z ~= 0
        fprintf('    Problema infactible \n');
        fprintf('Fi simplex primal \n'); 
        return;
    end
    
    fprintf('    Solució bàsica factible trobada, iteració  %2d \n', ite-1);
    
    % declaració de les variables de la faseII:

    for surt = find(VB > n)
        entra = VNB(1);
        [VNB, VB, invB] = actualitzar (entra, surt, VB, VNB, invB, -invB*A(:, entra));
    end
    
    VNB = VNB(VNB <= n); % treure les variables no-bàsiques majors que n
    
    % faseII del simplex:
    fprintf('  Fase II \n');
    [xb, tipus, VB, VNB, invB, ite] = fase(A, c, xb, VB, VNB, invB, ite, Bland);
    fprintf('  q = %2d,   rq = % 9.3f,   B(p) = %2d,   theta* = % 9.3f,   z = % 9.3f \n', 0, 0, 0, 0, c(VB)'*xb);
    
    if tipus == 1
        fprintf('    Problema factible il·limitat \n');
        fprintf('Fi simplex primal \n');
        return;
    end
    
    % i calcular la solució:
    z = c(VB)'*xb;
    
    fprintf('    Solució òptima trobada, iteració  %2d, z = % 9.6f \n', ite-1, z);
    fprintf('Fi simplex primal \n');
    
    fprintf('\n');
    fprintf('VB*= \n');
    disp(VB);
    fprintf('xb*= \n');
    disp(xb');
    fprintf('VNB*= \n');
    disp(VNB);
    fprintf('r*= \n');
    disp(c(VNB)' - c(VB)'*invB*A(:, VNB));
    fprintf('z*= \n');
    disp(z);
    
end


% la següent funció fa tant la faseI com la faseII
% aquí tipus serà 0 si ja s'ha trobat la solució òptima i 1 si és factible il·limitat
function [xb, tipus, VB, VNB, invB, ite] = fase(A, c, xb, VB, VNB, invB, ite, Bland)
    tipus = -1;
    while tipus == -1
        fprintf('    Iteració %2d : ', ite);
        ite = ite + 1;
        [xb, tipus, VB, VNB, invB] = iteracio(A, c, xb, VB, VNB, invB, Bland);
    end
end


% aquí el tipus serà -1 si s'ha fet una passa, 0 si ja s'ha trobat la solució òptima i 1 si és factible il·limitat
function [xb, tipus, VB, VNB, invB] = iteracio(A, c, xb, VB, VNB, invB, Bland)
    tipus = -1;

    r = c(VNB)' - c(VB)'*invB*A(:, VNB);
    r = r';
    if all(r >= -1e-8)
        tipus = 0;
        return;
    end
    
    % entra serà la variable que entra a la base
    if(Bland)
        entren = find(r < -1e-8);
    else
        entren = find(r == min(r));
    end
    entra = VNB(entren(1));
    
    % direcció de les variebles bàsiques
    db = -invB*A(:, entra);
    
    % comprovar si es il·limitat
    if all(db >= -1e-8)
        tipus = 1;
        return;
    end
    
    neg = db < -1e-8;
    theta = min(-xb(neg)./db(neg));
    
    % actualització de les variables:
    xb = xb + theta*db;
    
    % surt serà la posició de la variable que surt de la base
    surten = find(xb < 1e-8 & db < -1e-8); % possibles valors de surt
    if Bland
        surt = min(VB(surten));
        surt = find(VB == surt);
    else
        surt = surten(1);
    end
    
    xb(surt) = theta;
    
    % guardant els valors que s'hauran de mostrar
    q = entra;
    rq = r(VNB == entra);
    Bp = VB(surt);

    [VNB, VB, invB] = actualitzar (entra, surt, VB, VNB, invB, db);
    
    fprintf('  q = %2d,   rq = % 9.3f,   B(p) = %2d,   theta* = % 9.3f,   z = % 9.3f \n', q, rq, Bp, theta, c(VB)'*xb);

end



function [VNB, VB, invB] = actualitzar (entra, surt, VB, VNB, invB, db) % com aquestes variables s'havien s'actualizat en més d'un lloc del codi ho he fet com a funció
    VNB = VNB(VNB ~= entra);
    VNB = [VNB(VNB < VB(surt)), VB(surt), VNB(VNB > VB(surt))];  % mantenc VNB com a vector ordenat per a poder aplicar fàcilment la regla de Bland
    VB(surt) = entra;
    
    
    col_eta = -db/db(surt);
    fil_B = -invB(surt, :)/db(surt);
    invB = invB + col_eta*invB(surt, :);
    invB(surt, :) = fil_B;
end

