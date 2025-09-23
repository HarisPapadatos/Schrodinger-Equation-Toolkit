classdef Units
    properties (Constant)
        ev = 1;
        nm = 1;
        c  = 1;
        s  = 2.998e17;
        hbar = 197.326;
        me = 5.109e5;
        mm = 1000;
        um = 1e6;
        m  = 1e9;
        j  = 6.24e18;
        kg = 5.609e35;
        kev = 1000;
        mev = 1e6;
        evpc2 = 1;
        kevpc2 = 1000;
        mevpc2 = 1e6;
        hz = 3.336e-18;
    end

    methods (Static)
        function initialize()
            % persistent storage (class-wide, survives across calls)
            persistent keys_list values_list display_list valueMap displayMap

            if ~isempty(keys_list)
                warning('Units:AlreadyInitialized','Units already initialized.');
                return;
            end

            keys_list = {'ev','nm','c','s','hbar','me','mm','um','m','j','kg','kev','mev','evpc2','kevpc2','mevpc2','hz'};
            values_list = [ Units.ev, Units.nm, Units.c, Units.s, Units.hbar, Units.me, Units.mm, Units.um, Units.m, ...
                            Units.j, Units.kg, Units.kev, Units.mev, Units.evpc2, Units.kevpc2, Units.mevpc2, Units.hz ];
            display_list = {'eV','nm','c','s','\hbar','m_e','mm','\mu m','m','J','kg','keV','MeV','eV/c^2','keV/c^2','MeV/c^2','Hz'};
            valueMap   = containers.Map(keys_list, values_list);
            displayMap = containers.Map(keys_list, display_list);

            % stash into appdata (so other static methods can reach them)
            setappdata(0,'Units_keys',keys_list);
            setappdata(0,'Units_values',values_list);
            setappdata(0,'Units_display',display_list);
            setappdata(0,'Units_valueMap',valueMap);
            setappdata(0,'Units_displayMap',displayMap);
        end

        function v = value(key)
            if ~isappdata(0,'Units_valueMap')
                error('Units:NotInitialized','Call Units.initialize first.');
            end
            v = getappdata(0,'Units_valueMap');
            v = v(key);
        end

        function d = display(key)
            if ~isappdata(0,'Units_displayMap')
                error('Units:NotInitialized','Call Units.initialize first.');
            end
            d = getappdata(0,'Units_displayMap');
            d = d(key);
        end
    end
end
