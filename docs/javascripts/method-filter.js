const methodFilter = document.getElementById("method-filter");
if (methodFilter) {
    const input = document.getElementById("method-search");
    const status = document.getElementById("method-search-status");
    const groups = Array.from(methodFilter.parentElement.querySelectorAll("h2"), heading => {
        const content = [];
        for (let element = heading.nextElementSibling; element && element.tagName !== "H2"; element = element.nextElementSibling) {
            content.push(element);
        }
        return {heading, content, rows: content.flatMap(element => Array.from(element.querySelectorAll("tbody tr")))};
    }).filter(group => group.rows.length);

    function filterMethods() {
        const terms = input.value.toLowerCase().trim().split(/\s+/).filter(Boolean);
        let count = 0;
        for (const {heading, content, rows} of groups) {
            let matches = 0;
            for (const row of rows) {
                const text = heading.textContent.toLowerCase() + " " + row.textContent.toLowerCase();
                row.hidden = !terms.every(term => text.includes(term));
                if (!row.hidden) matches++;
            }
            heading.hidden = matches === 0;
            for (const element of content) element.hidden = matches === 0;
            count += matches;
        }
        status.textContent = terms.length ? `${count} matching entries` : "";
    }

    methodFilter.hidden = false;
    input.addEventListener("input", filterMethods);
}
