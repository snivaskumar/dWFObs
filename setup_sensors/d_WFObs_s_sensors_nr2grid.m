function [ d_grid,loc,element,d_strucObs ] = d_WFObs_s_sensors_nr2grid( stateid, Wp, d_Wp, d_strucObs)
%WFObs_s_sensors_nr2grid(stateid,Wp) Converts observation id to grid coordinates.

tur     = d_Wp{1}.tur;
loc     = cell(tur,1);
d_grid  = cell(tur,1);
if stateid <= 0
    error('State id too small. Check consistency with Wp.');
elseif stateid <= (Wp.Nx-3)*(Wp.Ny-2)
    element = 'u';
    grid.x = 3+floor( (stateid-1)/(Wp.Ny-2) );
    grid.y = 1+stateid-(grid.x-3)*(Wp.Ny-2);
    zz = [1:(Wp.Nx-3)*(Wp.Ny-2)];
    tmp1 = reshape(zz,(Wp.Ny-2),(Wp.Nx-3))';
    for i = 1:tur
        Nx = d_Wp{i}.mesh.Nx;
        Ny = d_Wp{i}.mesh.Ny;
        if ( ( d_Wp{i}.mesh.Nxb<=grid.x )&&( grid.x<=d_Wp{i}.mesh.Nxe ) )...
                &&( ( d_Wp{i}.mesh.Nyb<=grid.y )&&( grid.y<=d_Wp{i}.mesh.Nye ) )
            d_grid{i}.x = grid.x - d_Wp{i}.mesh.Nxb + 1;
            d_grid{i}.y = grid.y - d_Wp{i}.mesh.Nyb + 1;
            loc{i}.x = d_Wp{i}.mesh.ldxx2( d_grid{i}.x,1 );
            loc{i}.y = d_Wp{i}.mesh.ldyy(1,d_grid{i}.y);
            ZZ = [1:(Nx-3)*(Ny-2)];
            tmp = reshape(ZZ,(Ny-2),(Nx-3))';
            d_strucObs{i}.obs_array = tmp( d_grid{i}.x-3+1, d_grid{i}.y-2+1 );
%             d_strucObs{i}.obs_array_locu(end+1) = loc{i};
        end
    end
    Wp.ldxx2(grid.x,1);
    Wp.ldyy(1,grid.y);
    
elseif stateid <=(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)
    element = 'v';
    stateid  = stateid - (Wp.Nx-3)*(Wp.Ny-2);
    grid.x = 2+floor( (stateid-1)/(Wp.Ny-3));
    grid.y = 2+stateid-(grid.x-2)*(Wp.Ny-3);
    for i = 1:tur
        Nx = d_Wp{i}.mesh.Nx;
        Ny = d_Wp{i}.mesh.Ny;
        if ( ( d_Wp{i}.mesh.Nxb<=grid.x )&&( grid.x<=d_Wp{i}.mesh.Nxe ) )...
                &&( ( d_Wp{i}.mesh.Nyb<=grid.y )&&( grid.y<=d_Wp{i}.mesh.Nye ) )
            d_grid{i}.x = grid.x - d_Wp{i}.mesh.Nxb + 1;
            d_grid{i}.y = grid.y - d_Wp{i}.mesh.Nyb + 1;
            loc{i}.x = d_Wp{i}.mesh.ldxx( d_grid{i}.x,1 );
            loc{i}.y = d_Wp{i}.mesh.ldyy2(1,d_grid{i}.y);
            ZZ = [(Nx-3)*(Ny-2)+1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)];
            tmp = reshape(ZZ,(Ny-3),(Nx-2))';
            d_strucObs{i}.obs_array = tmp( d_grid{i}.x-3+1, d_grid{i}.y-2+1 );
%             d_strucObs{i}.obs_array_locv(end+1) = loc{i};
        end
    end
    Wp.ldxx(grid.x,1);
    Wp.ldyy2(1,grid.y);
    
elseif stateid <=(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2)-2
    element = 'p';
    stateid  = stateid - ((Wp.Nx-3)*(Wp.Ny-2) + (Wp.Nx-2)*(Wp.Ny-3));
    grid.x = 2+floor( (stateid-1)/(Wp.Ny-2));
    grid.y = 1+stateid-(grid.x-2)*(Wp.Ny-2);
    for i = 1:tur
        if ( ( d_Wp{i}.mesh.Nxb<grid.x )&&( grid.x<d_Wp{i}.mesh.Nxe ) )...
                &&( ( d_Wp{i}.mesh.Nyb<grid.y )&&( grid.y<d_Wp{i}.mesh.Nye ) )
            d_grid{i}.x = grid.x - d_Wp{i}.mesh.Nxb + 1;
            d_grid{i}.y = grid.y - d_Wp{i}.mesh.Nyb + 1;
            loc{i}.x = d_Wp{i}.mesh.ldxx2( d_grid{i}.x,1 );
            loc{i}.y = d_Wp{i}.mesh.ldyy(1,d_grid{i}.y);
            ZZ = [1:(d_Wp{i}.mesh.Nx-3)*(d_Wp{i}.mesh.Ny-2)];
            tmp = reshape(ZZ,(d_Wp{i}.mesh.Ny-2),(d_Wp{i}.mesh.Nx-3))';
            d_strucObs{i}.obs_array = tmp( d_grid{i}.x-3+1, d_grid{i}.y-2+1 );
        end
    end
    Wp.ldxx2(grid.x,1)
    Wp.ldyy(1,grid.y)
else
    error('State id too large. Check consistency with Wp.');
end;

end

