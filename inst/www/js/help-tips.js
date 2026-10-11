// Help icons (help_tip() in R/ui_modules.R): a Bootstrap tooltip on hover, on
// keyboard focus, and on a click or tap. A bare title attribute showed nothing
// on a click or a tap (user report, 2026-10-11).
//
// The icons sit inside input labels, and a click on a label is passed on to
// its input: a click on the Control group icon opened that menu instead. The
// click is stopped here, in the capture phase, before the label sees it.
$(function () {
  $('body').tooltip({
    selector: '.help-tip',
    container: 'body',
    placement: 'auto top',
    trigger: 'hover focus',
    template: '<div class="tooltip help-tip-tooltip" role="tooltip">' +
              '<div class="tooltip-arrow"></div><div class="tooltip-inner"></div></div>'
  });
  document.addEventListener('click', function (e) {
    var tip = e.target.closest ? e.target.closest('.help-tip') : null;
    if (!tip) return;
    e.preventDefault();
    e.stopPropagation();
    tip.focus();
  }, true);
});
