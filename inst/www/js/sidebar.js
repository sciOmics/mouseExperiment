/* Collapsible analysis sidebar panels
 *
 * Toggles the .sb-collapsed class on the nearest .collapsible-sidebar-row
 * ancestor. CSS handles hiding the .sidebar-panel-col and expanding
 * .main-panel-col to full width.
 */
$(document).ready(function () {
  $(document).on('click', '.sidebar-toggle-btn', function () {
    var $btn = $(this);
    var $row = $btn.closest('.collapsible-sidebar-row');

    if ($row.hasClass('sb-collapsed')) {
      $row.removeClass('sb-collapsed');
      $btn.html('<i class="fa fa-chevron-left fa-fw"></i> Collapse');
      $btn.attr('title', 'Collapse options panel');
    } else {
      $row.addClass('sb-collapsed');
      $btn.html('<i class="fa fa-chevron-right fa-fw"></i> Options');
      $btn.attr('title', 'Expand options panel');
    }
  });
});
