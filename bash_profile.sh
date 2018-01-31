# added by Anaconda2 2.5.0 installer
export PATH="/Users/Onur/anaconda/bin:$PATH"

# colors for ls command
export LS_COLORS='no=00:fi=00:di=01;44;37:ex=01;32:*.csv=00;36:*.tsv=00;36:*.tab=00;36:*.xml=00;94:*.xls=00;94:*.xlsx=00;94:*.json=00;94:*.hdf5=00;93:*.h5=00;93:*.pkl=00;93:*.pickle=00;93:*.ipynb=00;93:*.py=00;91:*.sh=00;91:*.slurm=00;42;30:*.sbatch=00;42;30:*.qsub=00;42;30:*.cmd=00;91:*.exe=01;91:*.com=01;91:*.bat=01;91:*.btm=01;91:*.dll=01;91:*.tar=00;31:*.tbz=00;31:*.tgz=00;31:*.rpm=00;31:*.deb=00;31:*.arj=00;31:*.taz=00;31:*.lzh=00;31:*.lzma=00;31:*.zip=00;31:*.zoo=00;31:*.z=00;31:*.Z=00;31:*.gz=00;31:*.bz2=00;31:*.tb2=00;31:*.tz2=00;31:*.tbz2=00;31:*.avi=01;35:*.bmp=01;35:*.fli=01;35:*.gif=01;35:*.jpg=01;35:*.jpeg=01;35:*.mng=01;35:*.mov=01;35:*.mpg=01;35:*.pcx=01;35:*.pbm=01;35:*.pgm=01;35:*.png=01;35:*.ppm=01;35:*.tga=01;35:*.tif=01;35:*.xbm=01;35:*.xpm=01;35:*.dl=01;35:*.gl=01;35:*.wmv=01;35:*.aiff=00;35:*.au=00;35:*.mid=00;35:*.mp3=00;35:*.ogg=00;35:*.voc=00;35:*.wav=00;35:*.svg=00;35:*.eps=00;35:*.pdf=00;35:'

# custom ls command
alias ll='ls -lrth --color'

# show git branch function
parse_git_branch() {
     git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'
}

# preferred git log
alias glog='git log --name-only --decorate --graph --all'

# nice extractor
extract () {
        if [ -f $1 ] ; then
          case $1 in
            *.tar.bz2)   tar xjf $1     ;;
            *.tar.gz)    tar xzf $1     ;;
            *.bz2)       bunzip2 $1     ;;
            *.rar)       unrar e $1     ;;
            *.gz)        gunzip $1      ;;
            *.tar)       tar xf $1      ;;
            *.tbz2)      tar xjf $1     ;;
            *.tgz)       tar xzf $1     ;;
            *.zip)       unzip $1       ;;
            *.Z)         uncompress $1  ;;
            *.7z)        7z x $1        ;;
            *)     echo "'$1' cannot be extracted via extract()" ;;
             esac
         else
             echo "'$1' is not a valid file"
         fi
}

# custom prompt
export PS1="\[\e[0;32m\][\u@\[\e[m\]\[\e[1;32m\]\h\[\e[m\]\[\e[0;32m\] ~ \@]\[\e[m\] \[\e[0;35m\]\w\[\e[m\]\[\e[0;33m\]\$(parse_git_branch)\[\e[m\] \[\e[0;32m\]>\[\e[m\] "

# ssh helpers
alias edison='ssh edison'
alias cori='ssh cori'
alias genepool='ssh genepool'
ssh-add ~/.ssh/id_rsa.nim
