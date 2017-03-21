https://ndownloader.figshare.com/articles/4640734/versions/1

for file in $(curl -s https://ndownloader.figshare.com/articles/4640734 | 
	grep href |
	sed 's/.*href="//' |
	sed 's/".*//' |
	grep '^[a-zA-Z].*'); do								curl -s -O https://ndownloader.figshare.com/articles/4640734/$file
	done
