BEGIN { status = 0; new_line_flag = false; synop_list = false; current_line = ""; } 

length($0) == 0 && status < 1000 {
    if (synop_list && length(current_line) > 0) {
        print current_line
    }
    status = 1000;
}

status == 50 {
    tmp = substr($0, index($0, $3))
    sub("<.*$", "", tmp)
    print "* " tmp
}

$2 == "Authors:" {
    if (synop_list) {
        if (length(current_line) > 0)
            print current_line
        synop_list = false
        print ""
    }
    status = 50
#    print ""
    print "== Authors =="
    print ""
}

substr($1, 0, 1) == "%" && (!new_line_flag || index($0, func_name) == 0) && status == 3 {
    if (NF == 1) {
        new_line_flag = 1
        synop_list = 1;
    } else {
        if (new_line_flag) {
            new_line_flag = false
            if (synop_list && length(current_line) > 0)
                print current_line
            current_line = "*"
        }
        text = substr($0, index($0, $2))
        gsub(func_name, tolower(func_name), text)
        if (synop_list) {
            current_line = current_line " " text
        } else {
            print text
        }
    }
}

substr($0, 0, 4) == "%   " && (status == 2 || status == 1) {
    status = 3
    new_line_flag = 1
    print ""
    print "== Synopsis =="
#    print ""
}

substr($0, 0, 1) == "%" && status == 1 {
    status = 2
    print ""
    print "== Description =="
    print ""
}

substr($1, 0, 1) == "%" && status == 2 {
    if (NR > 1) {
        if (substr($0, 2, 1) == " ") {
            print "=== " substr($0, index($0, $2)) " ==="
        } else {
            print substr($0, 2)
        }
    }
}

substr($1, 0, 1) == "%" && status == 0 {
    if (NR == 1) {
        func_name = toupper(substr(FILENAME, 0, length(FILENAME) - 2)); 
    }
    if (length($0) == 1) {
        status = 1;
    } else {
        print substr($0, 2, length($0) - 1)
    }
}

substr($1, 0, 1) == "%" && new_line_flag && index($0, func_name) > 0 && status == 3 {
    if (synop_list && length(current_line) > 0) {
        print current_line
    }
    new_line_flag = false
    line_start_ind = index($0, $2);
    func_ind = index($0, func_name);
    eq_ind = index($0, "=")
    header = "=== "
    if (eq_ind > 0 && eq_ind < func_ind) {
        tmp = substr($0, line_start_ind, eq_ind - line_start_ind)
        sub(" *$", "", tmp)
        header = header "''" tmp "'' = "
    }
    header = header tolower(func_name) "''" substr($0, index($0, "("), index($0, ")") - index($0, "(") + 1) "''"
    print ""
    print header " ==="
    desc_start = substr($0, index($0, ")") + 2)
    if (length(desc_start) > 0) {
        print ""
        sub(func_name, tolower(func_name), desc_start)
        print desc_start
    }
    synop_list = false;
    current_line = "";
}

{}
