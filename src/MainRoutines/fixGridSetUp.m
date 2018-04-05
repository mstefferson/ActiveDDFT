function systemObj = fixGridSetUp( systemObj, flags )
  % Change odd gridspacings to even unless it's one.
  if systemObj.n1 == 1
    systemObj.l1 = 1;
  else
    systemObj.n1 = systemObj.n1 + mod( systemObj.n1, 2 );
  end
  if systemObj.n2 == 1
    systemObj.l2 = 1;
  else
    systemObj.n2 = systemObj.n2 + mod( systemObj.n2, 2 );
  end
  if systemObj.n3 == 1
    systemObj.l3 = 1;
  else
    systemObj.n3 = systemObj.n3 + mod( systemObj.n3, 2 );
  end
  % Fix Ls if we want the box to be square
  if flags.SquareBox == 1
    systemObj.L_box = unique( [systemObj.l1 systemObj.l2] );
    systemObj.l1 = systemObj.L_box;
    systemObj.l2 = systemObj.L_box;
  end
  % Fix l1 is we want all Ns to be the same
  if flags.AllNsSame == 1
    if systemObj.n3 == 1
      Nvec = unique( [systemObj.n1 systemObj.n2] );
      systemObj.n1 = Nvec;  systemObj.n2 = Nvec;
    else
      Nvec = unique( [systemObj.n1 systemObj.n2 systemObj.n3] );
      systemObj.n1 = Nvec;  systemObj.n2 = Nvec;   systemObj.n3 = Nvec;
    end
  end
end

