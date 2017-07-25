# added by Anaconda2 2.5.0 installer
export PATH="/Users/Onur/anaconda/bin:$PATH"

# custom ls command
alias ll='ls -ltrhG'

# show git branch function
parse_git_branch() {
     git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'
}

# colors for ls command
export LSCOLORS='ExFxCxDxBxegedabagacad'

# custom prompt
export PS1="\[\e[0;32m\][\u@\[\e[m\]\[\e[1;32m\]\h\[\e[m\]\[\e[0;32m\] ~ \@]\[\e[m\] \[\e[0;35m\]\w\[\e[m\]\[\e[0;33m\]\$(parse_git_branch)\[\e[m\] \[\e[0;32m\]>\[\e[m\] "

# ssh helpers
alias edison='ssh edison'
alias cori='ssh cori'
alias genepool='ssh genepool'
ssh-add ~/.ssh/id_rsa.nim
