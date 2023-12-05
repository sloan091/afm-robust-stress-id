function  xp = getPrctile(x)

[~,srtidx] = sort(x);
dummy = 1:length(srtidx);
xp(srtidx) = 1 - dummy'/(length(srtidx) + 1);
xp = xp.';   
end