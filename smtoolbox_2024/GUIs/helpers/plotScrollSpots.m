function plotScrollSpots(spots)
f = figure('visible', 'off','position', [360 500 400 400]);
slhan = uicontrol('style', 'slider', 'position', [100 280 200 20],...
    'min',1,'max', 100,'callback', @callbackfn);
hsttext=uicontrol('style','text', 'position',[170 340 40 15],'visible','off');

axes('units','pixels','position',[100 50 200 200]);
movegui(f,'center')
set(f,'visible','on');

    function callbackfn(source,eventdata)
        i = get(slhan,'value');
        %set(hsttext,'visible','on','string',num2str(num))
        %imshow(spots(:,:,i));
        ax = gca;
        %ax.XLim=[0 2*pi]
    end
end