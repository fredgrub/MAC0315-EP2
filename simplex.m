function [ind v] = simplex(A,b,c,m,n,x,bind,Binv)
  optimal = false;
  count = 0;
  disp("Simplex: Fase 2\n")
  rule = -1;
  
  while !optimal
    # contador para as iterações
    printf("Iterando %i\n", count);
    count += 1;
  
    # imprime a solução viável básica atual
    printf("%i %f \n", [bind'; x(bind)']);
    printf("\nValor função objetivo: %f\n\n", c' * x);
    
    # índices das variáveis não-básicas
    indexes = [1:n];
    nonbasic = indexes(~ismember(indexes, bind));
    
    # vetor p usado para calcular os custos reduzidos
    p = (c(bind)' * Binv)';
    
    # custos reduzidos das variáveis não-basicas:
    # arrayfun mapeia a função que calcula o custo dado o índice
    # no vetor dos índices das variáveis não-básicas
    c_jbars = arrayfun(@(j) c(j) - (p' * A(:, j)), nonbasic);
    
    # imprime os custos reduzidos atuais
    printf("Custos reduzidos\n");
    printf("%i %f\n", [nonbasic; c_jbars]);
    
    # se nenhum custo reduzido é negativo, a solução ótima
    negative_c_jbars = c_jbars(c_jbars < 0);
    if isempty(negative_c_jbars)
      ind = 0;
      v = x;
      printf("\nSolução ótima encontrada com custo %f:\n", c' * x);
      printf("%i %f\n", [1:n; v']);
      return
    endif
    
    chosen_j = -1;    
    switch (rule)
      case 1 # custo reduzido mais negativo
        chosen_j = find(c == min(c))
      otherwise # regra de Bland
        for j = 1:(n-m)
          if c_jbars(j) < 0
            chosen_j = nonbasic(j);
            break;
          endif
        endfor
    endswitch
      
      
    endswitch
    
    # imprime o índice escolhido para entrar na base
    printf("\nEntra na base: %i\n\n", chosen_j);

    # vetor u
    u = Binv * A(:, chosen_j);
    
    # direção d
    d = zeros(n, 1);
    d(chosen_j) = 1;
    d(bind) = -u;
    
    # se u não tem componentes positivas, o custo ótimo é infinito na direção d
    if isempty(u(u > 0))
      ind = -1;
      v = d;
      printf("Custo ótimo é INFINITO na direção:\n");
      printf("%i %f\n", [1:n; v']);
      return
    endif
    
    # imprime as componentes básicas da direção básica
    printf("Direção\n");
    printf("%i %f\n", [bind'; d(bind)']);
    
    # agora procuramos theta*
    l = 0;
    thetastar = intmax;

    # varremos todos os índices com u(i) positivo e calculamos o x(bind(i))/u(i)
    # e selecionamos o menor
    # em caso de empate, escolhemos o i de menor bind(i)
    for i = 1:m
      if u(i) > 0
        theta = x(bind(i)) / u(i);
        if theta < thetastar || (abs(theta-thetastar) < eps && bind(i) < bind(l))
          thetastar = theta;
          l = i;
        endif
      endif
    endfor
    
    printf("\nTheta*\n");
    printf("%f\n\n", thetastar);
    
    # o indice l da base (que é o indice bind(l) da solucao) é escolhido pra sair
    printf("Sai da base: %i\n\n", bind(l));

    # nova solucao y
    y = zeros(n, 1);

    for i = 1:m
      if i == l
        continue
      endif;
      y(bind(i)) = x(bind(i)) - thetastar * u(i);
    endfor;

    y(chosen_j) = thetastar;

    # nova matriz inversa B
    tempB = [Binv u];

    # operacoes para arrumar a ultima coluna
    Q = eye(m);
    Q(:, l) = -u ./ u(l);
    Q(l, l) = 1 / u(l);

    newBinv = Q * tempB;

    # atualizamos as variaveis para a proxima iteracao
    Binv = newBinv(:, 1:m);
    x = y;
    bind(l) = chosen_j;
  endwhile
