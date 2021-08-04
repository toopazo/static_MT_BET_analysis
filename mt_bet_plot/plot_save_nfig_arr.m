function plot_save_nfig_arr(str1, str2, nfig_arr, savefig)      
    for nfig = nfig_arr       
        fig = figure(nfig);
        % set(fig,'units', 'centimeters', 'position', [0, 0, 20, 10]);
        
        filename = ['img/' str1 '_' str2 '_f' num2str(nfig) '.jpg'];
        saveas(fig, filename);     
        
        if savefig
            filename = ['img/' str1 '_' str2 '_f' num2str(nfig) '.fig'];
            saveas(fig, filename);
        end
    end
end

